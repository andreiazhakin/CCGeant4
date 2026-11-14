#include "G4CCElasticTableData.hh"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <mutex>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <cerrno>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace tabulated {
namespace {
[[noreturn]] void Fail(const std::string& s) { throw std::runtime_error("G4CCElasticTableData: "+s); }
using IsotopeKey = std::tuple<int,int,int,int>;
using ElementKey = std::pair<int,int>;
struct HeaderPrefix {
  char magic[8]; std::uint32_t version, header_bytes;
  std::int32_t zp,ap,zt,at; std::uint64_t energies,points;
  double target_mass,projectile_mass;
};
static_assert(sizeof(HeaderPrefix)==64,"Portable header prefix layout changed");
class Mapping final {
 public:
  explicit Mapping(const std::filesystem::path& p) {
#ifdef _WIN32
    file_=CreateFileW(p.c_str(),GENERIC_READ,FILE_SHARE_READ,nullptr,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,nullptr);
    if(file_==INVALID_HANDLE_VALUE) Fail("cannot open table "+p.string());
    LARGE_INTEGER n{}; if(!GetFileSizeEx(file_,&n)||n.QuadPart<256) { CloseHandle(file_); file_=INVALID_HANDLE_VALUE; Fail("invalid table size"); }
    bytes=static_cast<std::uint64_t>(n.QuadPart);
    mapping_=CreateFileMappingW(file_,nullptr,PAGE_READONLY,0,0,nullptr);
    if(mapping_) data=static_cast<const unsigned char*>(MapViewOfFile(mapping_,FILE_MAP_READ,0,0,0));
    if(!data) { if(mapping_) CloseHandle(mapping_); CloseHandle(file_); file_=INVALID_HANDLE_VALUE; mapping_=nullptr; Fail("cannot map "+p.string()); }
#else
    const int fd=::open(p.c_str(),O_RDONLY);
    if(fd<0) Fail("cannot open table "+p.string()+": "+std::strerror(errno));
    struct stat status{};
    if(::fstat(fd,&status)!=0 || status.st_size<256 ||
       static_cast<std::uintmax_t>(status.st_size)>std::numeric_limits<std::size_t>::max()) {
      ::close(fd); Fail("invalid table size");
    }
    bytes=static_cast<std::uint64_t>(status.st_size);
    void* mapped=::mmap(nullptr,static_cast<std::size_t>(bytes),PROT_READ,MAP_PRIVATE,fd,0);
    const int mapError=errno;
    ::close(fd);
    if(mapped==MAP_FAILED) Fail("cannot map "+p.string()+": "+std::strerror(mapError));
    data=static_cast<const unsigned char*>(mapped);
#endif
  }
  ~Mapping() {
#ifdef _WIN32
    if(data) UnmapViewOfFile(data);
    if(mapping_) CloseHandle(mapping_);
    if(file_!=INVALID_HANDLE_VALUE) CloseHandle(file_);
#else
    if(data) ::munmap(const_cast<unsigned char*>(data),static_cast<std::size_t>(bytes));
#endif
  }
  Mapping(const Mapping&)=delete; Mapping& operator=(const Mapping&)=delete;
  const unsigned char* data=nullptr; std::uint64_t bytes=0;
 private:
#ifdef _WIN32
  HANDLE file_=INVALID_HANDLE_VALUE,mapping_=nullptr;
#endif
};
double IntegralT(double t) {
  return std::abs(t)<1e-3 ? t*t*(.5+t*(-2./3.+t*(.75+t*(-.8+t*5./6.))))
                         : std::log1p(t)-t/(1.+t);
}
double Integral(double x0,double x1,double r0,double r1,double A) {
  const double d=x1-x0; if(d<=0.) return 0.;
  const double a=A+x0;
  return r0*d/(a*(a+d))+(r1-r0)/d*IntegralT(d/a);
}
void Probability(double p) { if(!std::isfinite(p)||p<0.||p>1.) Fail("probability outside [0,1]"); }
}

struct G4CCElasticTableData::Impl {
  std::vector<std::unique_ptr<Mapping>> maps;
  std::map<ElementKey,ProjectileTable> tables;
  std::set<ElementKey> indexed;
  std::uint64_t mapped_bytes=0;
};
G4CCElasticTableData::G4CCElasticTableData():impl_(new Impl) {}
G4CCElasticTableData::~G4CCElasticTableData()=default;
std::size_t G4CCElasticTableData::LoadedCount() const noexcept { return impl_->tables.size(); }
std::size_t G4CCElasticTableData::IndexedCount() const noexcept { return impl_->indexed.size(); }
std::uint64_t G4CCElasticTableData::MappedBytes() const noexcept { return impl_->mapped_bytes; }


std::shared_ptr<const G4CCElasticTableData>
G4CCElasticTableData::LoadFromEnvironment() {
  const char* configured = std::getenv("G4ELASTIC_TABLE_DIR");
  return LoadConfigured(configured ? configured : "");
}

std::shared_ptr<const G4CCElasticTableData>
G4CCElasticTableData::LoadConfigured(const std::string& configured) {
  if(configured.empty()) return {};
  std::filesystem::path supplied(configured);
  const auto directory = supplied.filename() == "tables.index.tsv" ? supplied.parent_path() : supplied;
  if(!std::filesystem::exists(directory/"tables.index.tsv")) return {};
  const auto path = std::filesystem::absolute(
      directory).lexically_normal().string();
  static std::mutex cacheMutex;
  static std::map<std::string, std::shared_ptr<const G4CCElasticTableData>> cache;
  std::lock_guard<std::mutex> guard(cacheMutex);
  auto found = cache.find(path);
  if (found != cache.end()) return found->second;
  auto data = Load(path);
  cache.emplace(path, data);
  return data;
}

std::shared_ptr<const G4CCElasticTableData> G4CCElasticTableData::Load(
    const std::string& directory,const std::vector<std::pair<int,int>>& selected) {
  const std::filesystem::path root(directory);
  const std::uint32_t endianProbe=1;
  if(*reinterpret_cast<const unsigned char*>(&endianProbe)!=1 || sizeof(double)!=8 ||
     !std::numeric_limits<double>::is_iec559)
    Fail("binary tables require little-endian IEEE-754 binary64 host");
  std::ifstream index(root/"tables.index.tsv"); if(!index) Fail("missing tables.index.tsv in "+directory);
  std::string line; std::getline(index,line);
  std::istringstream indexHeader(line);
  std::string marker, schema, extra; unsigned version = 0; std::size_t declaredCount = 0;
  if(!(indexHeader >> marker >> schema >> version >> declaredCount) ||
     marker != "#" || schema != "G4ELASTIC_TABLE_INDEX" || version != 1 ||
     declaredCount > 1000000 || (indexHeader >> extra)) Fail("unsupported table index/version/count");
  auto data=std::shared_ptr<G4CCElasticTableData>(new G4CCElasticTableData);
  while(std::getline(index,line)) {
    if(line.empty()) continue;
    std::istringstream row(line); std::string filename;
    int zp=0,ap=0,zt=0,at=0;
    if(!(row>>filename>>zp>>ap>>zt>>at)) Fail("malformed table index");
    if(filename.find_first_of("/\\:")!=std::string::npos) Fail("invalid index filename");
    const ElementKey key{zp,zt};
    if(!data->impl_->indexed.insert(key).second)
      Fail("ambiguous duplicate projectile-Z/target-Z table pair");
    const auto path=root/filename;
    if(!std::filesystem::is_regular_file(path)) Fail("indexed table missing: "+path.string());
    if(!selected.empty() && std::none_of(selected.begin(),selected.end(),
       [zp](const auto& isotope){return isotope.first==zp;})) continue;
    auto map=std::make_unique<Mapping>(path); HeaderPrefix h{}; std::memcpy(&h,map->data,sizeof(h));
    if(std::memcmp(h.magic,"G4ELT001",8)||h.version!=1||h.header_bytes!=256) Fail("unsupported binary header: "+path.string());
    if(IsotopeKey{h.zp,h.ap,h.zt,h.at}!=IsotopeKey{zp,ap,zt,at}) Fail("binary/index identity mismatch");
    if(h.energies<2||h.energies>1000000||h.points<2*h.energies||h.points>100000000) Fail("invalid energy/point counts");
    const std::uint64_t expected=256+8*(6*h.energies+1+5*h.points);
    if(map->bytes!=expected) Fail("truncated/extended binary table");
    if(!std::isfinite(h.target_mass)||!std::isfinite(h.projectile_mass)||h.target_mass<=0||h.projectile_mass<=0) Fail("invalid isotope mass");
    ProjectileTable t; t.projectile_Z=zp; t.projectile_A=ap; t.target_Z=zt; t.target_A=at;
    t.target_mass_MeV=h.target_mass; t.projectile_mass_MeV=h.projectile_mass;
    t.energy_count_=static_cast<std::size_t>(h.energies); t.point_count_=static_cast<std::size_t>(h.points);
    const unsigned char* cursor=map->data+256;
    auto doubles=[&cursor](std::size_t n) { const auto p=reinterpret_cast<const double*>(cursor); cursor+=8*n; return p; };
    t.energies_=doubles(t.energy_count_);
    t.offsets_=reinterpret_cast<const std::int64_t*>(cursor); cursor+=8*(t.energy_count_+1);
    t.sigma_=doubles(t.energy_count_); t.As_=doubles(t.energy_count_);
    t.p_cm_=doubles(t.energy_count_); t.C_=doubles(t.energy_count_);
    t.x_=doubles(t.point_count_); t.R_=doubles(t.point_count_); t.CDF_=doubles(t.point_count_);
    t.survival_=doubles(t.point_count_); t.cell_=doubles(t.point_count_);
    if(t.offsets_[0]!=0||t.offsets_[t.energy_count_]!=static_cast<std::int64_t>(t.point_count_)||!(t.energies_[0]>0.)) Fail("invalid table boundaries");
    for(std::size_t i=0;i<t.energy_count_;++i) {
      if(t.offsets_[i]<0 || t.offsets_[i+1]>static_cast<std::int64_t>(t.point_count_) ||
         t.offsets_[i+1]-t.offsets_[i]<2 || !std::isfinite(t.energies_[i]) ||
         (i && !(t.energies_[i]>t.energies_[i-1]))) Fail("invalid offsets/energy grid");
      for(double v:{t.sigma_[i],t.As_[i],t.p_cm_[i],t.C_[i]}) if(!std::isfinite(v)||v<=0.) Fail("invalid node scalar");
      const auto n=t.Node(i);
      if(n.x[0]!=0.||n.x[n.size-1]!=1.||n.CDF[0]!=0.||n.CDF[n.size-1]!=1.||n.survival[0]!=1.||n.survival[n.size-1]!=0.||n.cell_mb[n.size-1]!=0.) Fail("invalid node endpoints");
    }
    data->impl_->mapped_bytes+=map->bytes;
    data->impl_->maps.push_back(std::move(map)); data->impl_->tables.emplace(key,t);
  }
  if(data->impl_->indexed.size()!=declaredCount) Fail("incomplete table index");
  return data;
}
const ProjectileTable* G4CCElasticTableData::FindByZ(int zp,int zt) const {
  const ElementKey key{zp,zt}; const auto i=impl_->tables.find(key);
  if(i!=impl_->tables.end()) return &i->second;
  if(impl_->indexed.count(key)) Fail("covered projectile was not preloaded; load all tables for hadronic residues");
  return nullptr;
}
const ProjectileTable* G4CCElasticTableData::Find(int zp,int,int zt,int) const {
  return FindByZ(zp,zt);
}
NodeView ProjectileTable::Node(std::size_t i) const {
  if(i>=energy_count_) Fail("node index out of bounds");
  const auto start=static_cast<std::size_t>(offsets_[i]);
  NodeView n; n.size=static_cast<std::size_t>(offsets_[i+1]-offsets_[i]);
  n.energy_MeVu=energies_[i]; n.sigma_mb=sigma_[i]; n.As=As_[i]; n.C_mb=C_[i]; n.p_cm_MeV=p_cm_[i];
  n.x=x_+start; n.R=R_+start; n.CDF=CDF_+start; n.survival=survival_+start; n.cell_mb=cell_+start; return n;
}
EnergyBracket ProjectileTable::Lookup(double e) const {
  EnergyBracket b;
  if(!std::isfinite(e)||e<energies_[0]||e>energies_[energy_count_-1]) return b;
  const auto p=std::lower_bound(energies_,energies_+energy_count_,e);
  b.valid=true; b.hi=static_cast<std::size_t>(p-energies_);
  if(*p==e) { b.lo=b.hi; b.sigma_mb=sigma_[b.lo]; return b; }
  b.lo=b.hi-1; b.w=(e-energies_[b.lo])/(energies_[b.hi]-energies_[b.lo]);
  b.sigma_mb=(1.-b.w)*sigma_[b.lo]+b.w*sigma_[b.hi]; return b;
}
std::size_t ProjectileTable::SelectNode(const EnergyBracket& b,double u) const {
  Probability(u); if(!b.valid) Fail("cannot sample an uncovered energy");
  return b.lo==b.hi || u>=b.w*sigma_[b.hi]/b.sigma_mb ? b.lo:b.hi;
}
double ProjectileTable::SampleX(const EnergyBracket& b,double u_node,double u) const { return SampleNodeX(SelectNode(b,u_node),u); }
double ProjectileTable::SampleNodeX(std::size_t i,double u) const {
  Probability(u); return u>.5 ? InverseNode(i,1.-u,true):InverseNode(i,u,false);
}
double ProjectileTable::InverseNode(std::size_t i,double u,bool tail) const {
  Probability(u); if(u==0.) return tail?1.:0.; if(u==1.) return tail?0.:1.;
  const auto n=Node(i); std::size_t j=0;
  if(tail) j=static_cast<std::size_t>(std::upper_bound(n.survival,n.survival+n.size,u,std::greater<double>())-n.survival-1);
  else j=static_cast<std::size_t>(std::upper_bound(n.CDF,n.CDF+n.size,u)-n.CDF-1);
  j=std::min(j,n.size-2);
  const double x0=n.x[j],d=n.x[j+1]-x0,a=n.As+x0,r0=n.R[j],slope=(n.R[j+1]-r0)/d;
  const double cell=n.cell_mb[j]/n.C_mb;
  const double remaining=tail?std::max(0.,u-n.survival[j+1])*n.sigma_mb/n.C_mb:0.;
  const double target=std::clamp(tail?cell-remaining:std::max(0.,u-n.CDF[j])*n.sigma_mb/n.C_mb,0.,cell);
  if(target==0.) return x0; if(target==cell) return x0+d;
  const double fraction=target/cell;
  double lo=0.,hi=d,delta=a*d*fraction/(a+d*(1.-fraction));
  const bool reverse=tail && remaining<.5*cell;
  for(int k=0;k<64;++k) {
    const double r=r0+slope*delta;
    const double left=r0*delta/(a*(a+delta))+slope*IntegralT(delta/a);
    const double residual=reverse?remaining-Integral(x0+delta,x0+d,r,n.R[j+1],n.As):left-target;
    const double scale=reverse?remaining:target;
    if(std::abs(residual)<=2e-14*std::max(scale,cell*1e-14)) return x0+delta;
    if(residual<0.) lo=delta; else hi=delta;
    if(hi-lo<=4.*std::numeric_limits<double>::epsilon()*std::max(std::abs(delta),d*1e-6)) return x0+.5*(lo+hi);
    const double derivative=r/((a+delta)*(a+delta));
    const double proposed=delta-residual/derivative;
    delta=(k<12 && std::isfinite(proposed) && proposed>lo && proposed<hi)?proposed:.5*(lo+hi);
  }
  return x0+.5*(lo+hi);
}
double ProjectileTable::CumulativeNode(std::size_t i,double x,bool tail) const {
  Probability(x); if(x==0.) return tail?1.:0.; if(x==1.) return tail?0.:1.;
  const auto n=Node(i);
  const auto j=std::min(static_cast<std::size_t>(std::upper_bound(n.x,n.x+n.size,x)-n.x-1),n.size-2);
  const double r=n.R[j]+(n.R[j+1]-n.R[j])*(x-n.x[j])/(n.x[j+1]-n.x[j]);
  return tail?n.survival[j+1]+n.C_mb*Integral(x,n.x[j+1],r,n.R[j+1],n.As)/n.sigma_mb
             :n.CDF[j]+n.C_mb*Integral(n.x[j],x,n.R[j],r,n.As)/n.sigma_mb;
}
double ProjectileTable::Cumulative(const EnergyBracket& b,double x,bool tail) const {
  if(!b.valid) Fail("CDF requested outside table coverage");
  if(b.lo==b.hi) return CumulativeNode(b.lo,x,tail);
  return ((1.-b.w)*sigma_[b.lo]*CumulativeNode(b.lo,x,tail)+b.w*sigma_[b.hi]*CumulativeNode(b.hi,x,tail))/b.sigma_mb;
}
}
