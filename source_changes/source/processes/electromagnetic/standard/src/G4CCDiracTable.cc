#include "G4CCDiracTable.hh"
#include "G4Exception.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double kT0 = G4CCDiracTable::TabulatedMinTransfer();
constexpr double kTMin = G4CCDiracTable::MinTransfer();
constexpr double kElectronMassMeV = 0.51099895;
constexpr double kAtomicMassMeV = 931.494;
using Moments = std::array<double, 3>;

[[noreturn]] void Fail(const std::string& text)
{
  G4Exception("G4CCDiracTable", "CCDiracTable001", FatalException, text.c_str());
  throw std::runtime_error(text);
}

bool Close(double a, double b)
{
  return std::abs(a-b) <= 1.e-10*std::max(std::abs(a), std::abs(b));
}

std::vector<std::filesystem::path> Files(const G4String& input)
{
  if(input.empty()) { Fail("Required Dirac CSV path was not supplied."); }
  const std::filesystem::path path(input.c_str());
  std::vector<std::filesystem::path> files;
  if(std::filesystem::is_directory(path)) {
    for(const auto& entry : std::filesystem::directory_iterator(path)) {
      if(entry.is_regular_file() && entry.path().extension() == ".csv") {
        files.push_back(entry.path());
      }
    }
    std::sort(files.begin(), files.end());
  } else if(std::filesystem::is_regular_file(path)) {
    files.push_back(path);
  }
  if(files.empty()) { Fail("No CSV files found at: " + std::string(input)); }
  return files;
}

void ReadHeader(std::ifstream& stream, const std::string& expected,
                const std::filesystem::path& path)
{
  std::string line;
  if(!std::getline(stream, line)) { Fail("Missing CSV header: " + path.string()); }
  if(!line.empty() && line.back() == '\r') { line.pop_back(); }
  if(line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xef &&
     static_cast<unsigned char>(line[1]) == 0xbb &&
     static_cast<unsigned char>(line[2]) == 0xbf) { line.erase(0,3); }
  if(line != expected) { Fail("Unsupported CSV header: " + path.string()); }
}

template<std::size_t N> std::array<double, N> Parse(const std::string& line)
{
  std::array<double, N> result{};
  const char* cursor = line.c_str();
  for(std::size_t i=0; i<N; ++i) {
    char* end = nullptr;
    errno = 0;
    result[i] = std::strtod(cursor, &end);
    if(end == cursor || errno == ERANGE || !std::isfinite(result[i])) {
      Fail("Invalid/nonfinite numeric value in Dirac CSV.");
    }
    if(i+1 < N) {
      if(*end != ',') { Fail("Wrong Dirac CSV column count."); }
      cursor = end+1;
    } else {
      while(*end == ' ' || *end == '\t' || *end == '\r') { ++end; }
      if(*end != '\0') { Fail("Unexpected text after Dirac CSV row."); }
    }
  }
  return result;
}

int Integer(double x, int low, int high)
{
  if(x < low || x > high || std::floor(x) != x) { Fail("Invalid integer in Dirac CSV."); }
  return static_cast<int>(x);
}

double CellMoment(double t, double d, double upper, double slope, int order)
{
  if(upper <= t) { return 0.; }
  const double l = std::log1p((upper-t)/t);
  const double v = (slope + order + 1.)*l;
  const double exprel = std::abs(v) < 1.e-12 ? 1.+0.5*v+v*v/6. : std::expm1(v)/v;
  return d*std::pow(t, order+1)*l*exprel;
}

struct Spectrum {
  std::vector<double> density;
  std::array<std::vector<double>, 3> tail;
};

struct EnergyNode {
  double energy=0.;
  std::vector<double> t;
  std::vector<Spectrum> spectra;
};

struct Projectile {
  int nmax=-1;
  std::vector<EnergyNode> nodes;
};

struct StopNode {
  double energy=0., b=0., bLeft=0., atima=0., atimaLeft=0.;
};

struct Bracket {
  const EnergyNode* lo=nullptr;
  const EnergyNode* hi=nullptr;
  double weight=0.;
};

Bracket FindBracket(const Projectile& projectile, double energy)
{
  const auto& nodes = projectile.nodes;
  if(!std::isfinite(energy) || energy < nodes.front().energy || energy > nodes.back().energy) {
    Fail("Requested Dirac energy lies outside supplied CSV coverage.");
  }
  const auto hi = std::lower_bound(nodes.begin(), nodes.end(), energy,
    [](const EnergyNode& n, double e) { return n.energy < e; });
  if(hi->energy == energy) { return {&*hi, &*hi, 0.}; }
  const auto lo = hi-1;
  return {&*lo, &*hi, std::log(energy/lo->energy)/std::log(hi->energy/lo->energy)};
}

Moments NodeTail(const EnergyNode& node, int n, double cut)
{
  const auto& spectrum = node.spectra[static_cast<std::size_t>(n)];
  Moments result{};
  if(cut >= node.t.back()) { return result; }
  if(cut <= node.t.front()) {
    for(int k=0; k<3; ++k) { result[k] = spectrum.tail[k].front(); }
    return result;
  }
  const auto upper = std::upper_bound(node.t.begin(), node.t.end(), cut);
  const std::size_t i = static_cast<std::size_t>(upper-node.t.begin()-1);
  const double slope = std::log(spectrum.density[i+1]/spectrum.density[i]) /
                       std::log(node.t[i+1]/node.t[i]);
  const double density = spectrum.density[i]*std::exp(slope*std::log(cut/node.t[i]));
  for(int k=0; k<3; ++k) {
    result[k] = spectrum.tail[k][i+1] + CellMoment(cut, density, node.t[i+1], slope, k);
  }
  return result;
}

double NodeDensity(const EnergyNode& node, int n, double t)
{
  if(t < node.t.front() || t > node.t.back()) { return 0.; }
  const auto& d = node.spectra[static_cast<std::size_t>(n)].density;
  if(t == node.t.back()) { return d.back(); }
  const auto upper = std::upper_bound(node.t.begin(), node.t.end(), t);
  const std::size_t i = static_cast<std::size_t>(upper-node.t.begin()-1);
  const double w = std::log(t/node.t[i])/std::log(node.t[i+1]/node.t[i]);
  return std::exp((1.-w)*std::log(d[i])+w*std::log(d[i+1]));
}

double NodeSample(const EnergyNode& node, int n, double cut, double u)
{
  if(u <= 0.) { return std::max(cut, node.t.front()); }
  if(u >= 1.) { return node.t.back(); }
  const auto& spectrum = node.spectra[static_cast<std::size_t>(n)];
  const auto& tail = spectrum.tail[0];
  const double survival = (1.-u)*NodeTail(node,n,cut)[0];
  const auto next = std::upper_bound(tail.begin(),tail.end(),survival,std::greater<double>());
  const std::size_t i = std::min(static_cast<std::size_t>(next-tail.begin()-1),node.t.size()-2);
  const double fraction = std::clamp((tail[i]-survival)/(tail[i]-tail[i+1]),0.,1.);
  if(fraction == 0.) { return std::max(cut,node.t[i]); }
  if(fraction == 1.) { return node.t[i+1]; }
  const double l = std::log(node.t[i+1]/node.t[i]);
  const double power = std::log(spectrum.density[i+1]/spectrum.density[i])/l + 1.;
  const double logRatio = std::abs(power*l)<1.e-10 ? fraction*l :
    std::log1p(fraction*std::expm1(power*l))/power;
  return std::clamp(node.t[i]*std::exp(logRatio),std::max(cut,node.t[i]),node.t[i+1]);
}

struct MappedNode {
  const EnergyNode& node;
  int n;
  double a, b, top;

  MappedNode(const EnergyNode& source, int electrons, double maximum)
    : node(source), n(electrons),
      a((maximum-kT0)/(source.t.back()-kT0)), b(kT0*(1.-a)), top(maximum) {}

  double FirstSlope() const {
    const auto& d=node.spectra[static_cast<std::size_t>(n)].density;
    return std::log(d[1]/d[0]) / std::log1p(a*(node.t[1]-kT0)/kT0);
  }

  double FirstDensity() const {
    return node.spectra[static_cast<std::size_t>(n)].density[0]/a;
  }

  double NativeCut(double cut) const {
    return std::clamp(kT0+(std::max(cut,kT0)-kT0)/a,kT0,node.t.back());
  }

  Moments Extension(double cut) const {
    Moments result{};
    if(cut>=kT0) {return result;}
    const double slope=FirstSlope();
    const double d=FirstDensity()*std::exp(slope*std::log(cut/kT0));
    for(int k=0;k<3;++k) {result[k]=CellMoment(cut,d,kT0,slope,k);}
    return result;
  }

  Moments Tail(double cut) const {
    if(cut>=top) {return {};}
    const auto v=NodeTail(node,n,NativeCut(cut));
    auto result=Extension(cut);
    result[0]+=v[0];
    result[1]+=a*v[1]+b*v[0];
    result[2]+=a*a*v[2]+2.*a*b*v[1]+b*b*v[0];
    return result;
  }

  double Density(double t) const {
    if(t<kTMin||t>top) {return 0.;}
    if(t<kT0) {return FirstDensity()*std::exp(FirstSlope()*std::log(t/kT0));}
    return NodeDensity(node,n,NativeCut(t))/a;
  }

  double Sample(double cut,double u) const {
    if(u<=0.) {return cut;}
    if(u>=1.) {return top;}
    const double low=Extension(cut)[0];
    const double nativeCut=NativeCut(cut);
    const double high=NodeTail(node,n,nativeCut)[0];
    const double target=u*(low+high);
    if(low>0. && target<low) {
      const double fraction=std::clamp(target/low,0.,1.);
      const double l=std::log(kT0/cut), power=FirstSlope()+1.;
      const double logRatio=std::abs(power*l)<1.e-10 ? fraction*l :
        std::log1p(fraction*std::expm1(power*l))/power;
      return std::clamp(cut*std::exp(logRatio),cut,kT0);
    }
    const double nativeU=low>0. ? std::clamp((target-low)/high,0.,1.) : u;
    return std::clamp(kT0+a*(NodeSample(node,n,nativeCut,nativeU)-kT0),cut,top);
  }
};

double QueryMaximum(double energy,double maximum) {
  const double tabulated=G4CCDiracTable::MaxTransfer(energy);
  if(maximum==0.) {return tabulated;}
  if(!std::isfinite(maximum)||maximum<=kT0||maximum>tabulated*(1.+1.e-12)) {
    Fail("Invalid maximum transfer for mapped Dirac spectrum.");
  }
  return maximum;
}

void CheckQuery(int n, double cut)
{
  if(n < 0 || n > 27) { Fail("Dirac tables support electron counts N=0..27."); }
  if(!std::isfinite(cut) || cut < kTMin) {
    Fail("Dirac electron production threshold must be at least 0.99 keV.");
  }
}

double ScaleForN(int z, int n, int nmax)
{
  if(n <= nmax) { return 1.; }
  const double ratio = static_cast<double>(z-n)/static_cast<double>(z-nmax);
  return ratio*ratio;
}
}

struct G4CCDiracTable::Impl {
  std::map<int, Projectile> projectiles;
  std::map<std::pair<int,int>, std::vector<StopNode>> stopping;
  double minimum=30., maximum=2000.;
  std::size_t memory=0;

  const Projectile& Get(int z) const {
    const auto found=projectiles.find(z);
    if(found==projectiles.end()) { Fail("No Dirac projectile CSV for Z="+std::to_string(z)); }
    return found->second;
  }

  const std::vector<StopNode>& GetStop(int z,int zt) const {
    const auto found=stopping.find({z,zt});
    if(found==stopping.end()) {
      Fail("No Dirac stopping CSV curve for Z="+std::to_string(z)+", Zt="+std::to_string(zt));
    }
    return found->second;
  }

  double InterpolateStop(int z,int zt,double energy,bool atima,bool left) const {
    const auto& curve=GetStop(z,zt);
    if(!std::isfinite(energy)||energy<curve.front().energy||energy>curve.back().energy) {
      Fail("Requested stopping energy lies outside supplied CSV coverage.");
    }
    const auto hi=std::lower_bound(curve.begin(),curve.end(),energy,
      [](const StopNode& p,double e){return p.energy<e;});
    if(hi->energy==energy) {
      return atima ? (left?hi->atimaLeft:hi->atima) : (left?hi->bLeft:hi->b);
    }
    const auto lo=hi-1;
    const double w=std::log1p((energy-lo->energy)/lo->energy)/
                   std::log1p((hi->energy-lo->energy)/lo->energy);
    return atima ? (1.-w)*lo->atima+w*hi->atimaLeft : (1.-w)*lo->b+w*hi->bLeft;
  }

  void LoadDelta(const G4String& path) {
    std::string header="Zp,Nmax,E_MeV_u,T_MeV";
    for(int n=0;n<=15;++n) {header+=",D"+std::to_string(n)+"_cm2_MeV";}
    for(const auto& file:Files(path)) {
      std::ifstream stream(file);
      if(!stream) {Fail("Cannot open Dirac CSV: "+file.string());}
      ReadHeader(stream,header,file);
      std::string line;
      std::size_t row=1;
      while(std::getline(stream,line)) {
        ++row;
        if(line.empty()||line=="\r") {continue;}
        try {
          const auto values=Parse<20>(line);
          const int z=Integer(values[0],29,92),nmax=Integer(values[1],0,15);
          if(nmax!=(z<=50?10:15)) {Fail("Unexpected calculated Dirac Nmax.");}
          auto& projectile=projectiles[z];
          if(projectile.nmax<0) {projectile.nmax=nmax;}
          if(projectile.nmax!=nmax) {Fail("Inconsistent Dirac Nmax for one projectile.");}
          if(projectile.nodes.empty()||projectile.nodes.back().energy!=values[2]) {
            if(!projectile.nodes.empty()&&values[2]<=projectile.nodes.back().energy) {
              Fail("Dirac energy blocks must be strictly increasing for each Z.");
            }
            EnergyNode node;
            node.energy=values[2];
            node.spectra.resize(static_cast<std::size_t>(nmax+1));
            projectile.nodes.push_back(std::move(node));
          }
          auto& node=projectile.nodes.back();
          if(values[3]<kT0||(!node.t.empty()&&values[3]<=node.t.back())) {
            Fail("Dirac transfer grid must be positive and strictly increasing.");
          }
          node.t.push_back(values[3]);
          for(int n=0;n<=15;++n) {
            if(n<=nmax) {
              if(!(values[4+n]>0.)) {Fail("Dirac differential cross sections must be positive.");}
              node.spectra[static_cast<std::size_t>(n)].density.push_back(values[4+n]);
            } else if(values[4+n]!=0.) {Fail("Unused Dirac CSV columns must be zero.");}
          }
        } catch(const std::exception& e) {
          Fail(file.string()+":"+std::to_string(row)+": "+e.what());
        }
      }
      if(stream.bad()) {Fail("I/O failure reading "+file.string());}
    }
    if(projectiles.empty()) {Fail("No Dirac projectile spectra loaded.");}
    for(auto& entry:projectiles) {
      auto& projectile=entry.second;
      if(projectile.nodes.size()<2||projectile.nodes.front().energy!=minimum||
         projectile.nodes.back().energy!=maximum) {Fail("Dirac CSV must cover 30..2000 MeV/u.");}
      for(auto& node:projectile.nodes) {
        if(node.t.size()<2||node.t.front()!=kT0||
           !Close(node.t.back(),G4CCDiracTable::MaxTransfer(node.energy))) {
          Fail("Dirac transfer endpoints do not match the documented kinematics.");
        }
        node.t.shrink_to_fit();
        memory+=node.t.capacity()*sizeof(double);
        for(auto& spectrum:node.spectra) {
          spectrum.density.shrink_to_fit();
          memory+=spectrum.density.capacity()*sizeof(double);
          for(auto& tail:spectrum.tail) {tail.resize(node.t.size(),0.);memory+=tail.capacity()*sizeof(double);}
          for(std::size_t i=node.t.size()-1;i-->0;) {
            const double slope=std::log(spectrum.density[i+1]/spectrum.density[i])/
                               std::log(node.t[i+1]/node.t[i]);
            for(int k=0;k<3;++k) {
              const double cell=CellMoment(node.t[i],spectrum.density[i],node.t[i+1],slope,k);
              spectrum.tail[k][i]=spectrum.tail[k][i+1]+cell;
              if(!(cell>0.)||!std::isfinite(spectrum.tail[k][i])) {Fail("Invalid Dirac moment integral.");}
            }
          }
        }
      }
    }
  }

  void LoadStopping(const G4String& path) {
    for(const auto& file:Files(path)) {
      std::ifstream stream(file);
      if(!stream) {Fail("Cannot open stopping CSV: "+file.string());}
      ReadHeader(stream,"Zp,Zt,E_MeV_u,B_MeV_cm2,B_left_MeV_cm2,ATIMA_MeV_cm2,ATIMA_left_MeV_cm2",file);
      std::string line;
      std::size_t row=1;
      while(std::getline(stream,line)) {
        ++row;
        if(line.empty()||line=="\r") {continue;}
        try {
          const auto values=Parse<7>(line);
          const int z=Integer(values[0],29,92),zt=Integer(values[1],1,96);
          if(zt>=85&&zt<=87) {Fail("Unsupported target Z in stopping CSV.");}
          if(projectiles.find(z)==projectiles.end()) {Fail("Stopping projectile has no matching Dirac CSV.");}
          auto& curve=stopping[{z,zt}];
          if(!curve.empty()&&values[2]<=curve.back().energy) {Fail("Stopping energies must increase strictly.");}
          if(values[3]<0.||values[4]<0.||!(values[5]>0.)||!(values[6]>0.)) {
            Fail("Invalid soft coefficient or equilibrium stopping.");
          }
          curve.push_back({values[2],values[3],values[4],values[5],values[6]});
        } catch(const std::exception& e) {
          Fail(file.string()+":"+std::to_string(row)+": "+e.what());
        }
      }
      if(stream.bad()) {Fail("I/O failure reading "+file.string());}
    }
    for(auto& entry:stopping) {
      auto& curve=entry.second;
      if(curve.size()<2||curve.front().energy!=minimum||curve.back().energy!=maximum) {
        Fail("Stopping CSV must cover 30..2000 MeV/u for every supplied pair.");
      }
      curve.shrink_to_fit();memory+=curve.capacity()*sizeof(StopNode);
    }
    for(const auto& projectile:projectiles) {
      const auto first=stopping.lower_bound({projectile.first,0});
      if(first==stopping.end()||first->first.first!=projectile.first) {
        Fail("Dirac projectile has no supplied elemental stopping curve: Z="+
             std::to_string(projectile.first));
      }
    }
  }
};

G4CCDiracTable::G4CCDiracTable()=default;
G4CCDiracTable::~G4CCDiracTable()=default;
G4CCDiracTable& G4CCDiracTable::Instance()
{
  static G4CCDiracTable table;
  return table;
}

void G4CCDiracTable::Configure(const G4String& deltaPath,const G4String& stoppingPath)
{
  std::lock_guard<std::mutex> guard(fMutex);
  if(deltaPath==fDeltaPath&&stoppingPath==fStoppingPath) {return;}
  if(fLoaded.load(std::memory_order_acquire)||G4Threading::IsWorkerThread()||
     G4StateManager::GetStateManager()->GetCurrentState()!=G4State_PreInit) {
    Fail("Dirac CSV paths can only be configured on the master before initialization.");
  }
  if(deltaPath.empty()||stoppingPath.empty()) {Fail("Both delta and stopping CSV paths are required.");}
  fDeltaPath=deltaPath;fStoppingPath=stoppingPath;
}

void G4CCDiracTable::Validate() const
{
  if(fLoaded.load(std::memory_order_acquire)) {return;}
  std::lock_guard<std::mutex> guard(fMutex);
  if(fLoaded.load(std::memory_order_acquire)) {return;}
  auto data=std::make_unique<Impl>();
  data->LoadDelta(fDeltaPath);
  data->LoadStopping(fStoppingPath);
  fData=std::move(data);
  fLoaded.store(true,std::memory_order_release);
}

G4bool G4CCDiracTable::HasProjectile(G4int z) const
{
  Validate();return fData->projectiles.find(z)!=fData->projectiles.end();
}
G4bool G4CCDiracTable::HasPair(G4int z,G4int zt) const
{
  Validate();return fData->stopping.find({z,zt})!=fData->stopping.end();
}
G4double G4CCDiracTable::EnergyMin() const {Validate();return fData->minimum;}
G4double G4CCDiracTable::EnergyMax() const {Validate();return fData->maximum;}
std::size_t G4CCDiracTable::MemoryBytes() const {Validate();return fData->memory;}
G4double G4CCDiracTable::MaxTransfer(G4double energy)
{
  const double tau=energy/kAtomicMassMeV;
  return 2.*kElectronMassMeV*tau*(tau+2.);
}

std::array<G4double,3> G4CCDiracTable::DeltaMoments(G4int z,G4int n,G4double energy,
                                                G4double cut,G4double maximum) const
{
  CheckQuery(n,cut);Validate();
  const auto& projectile=fData->Get(z);
  const auto bracket=FindBracket(projectile,energy);
  Moments result{};
  const double top=QueryMaximum(energy,maximum);
  if(cut>=top) {return result;}
  const int nn=std::min(n,projectile.nmax);
  for(int j=0;j<2;++j) {
    const auto& node=*(j==0?bracket.lo:bracket.hi);
    const double weight=j==0?1.-bracket.weight:bracket.weight;
    if(weight==0.) {continue;}
    const auto tail=MappedNode(node,nn,top).Tail(cut);
    for(int k=0;k<3;++k) {result[k]+=weight*tail[k];}
  }
  const double scale=ScaleForN(z,n,projectile.nmax);
  for(auto& value:result) {value*=scale;}
  return result;
}

G4double G4CCDiracTable::DeltaCrossSectionUpperBound(G4int z,G4int n,
  G4double energyLow,G4double energyHigh,G4double cut,
  G4double maximumLow,G4double maximumHigh) const
{
  CheckQuery(n,cut);Validate();
  const auto& projectile=fData->Get(z);
  const auto low=FindBracket(projectile,energyLow);
  const auto high=FindBracket(projectile,energyHigh);
  const double topLow=QueryMaximum(energyLow,maximumLow);
  const double topHigh=QueryMaximum(energyHigh,maximumHigh);
  if(cut>=topHigh) {return 0.;}
  const int nn=std::min(n,projectile.nmax);
  double bound=0.;
  for(const auto* node=low.lo;node<=high.hi;++node) {
    double value=MappedNode(*node,nn,topHigh).Tail(std::max(cut,kT0))[0];
    if(cut<kT0) {
      const MappedNode mapped(*node,nn,topLow);
      const double slope=std::min(0.,mapped.FirstSlope());
      value+=(kT0-cut)*mapped.FirstDensity()*std::exp(slope*std::log(cut/kT0));
    }
    bound=std::max(bound,value);
  }
  return bound*ScaleForN(z,n,projectile.nmax)*(1.+1.e-12);
}

G4double G4CCDiracTable::DeltaDensity(G4int z,G4int n,G4double energy,G4double t,
                                     G4double maximum) const
{
  CheckQuery(n,kT0);Validate();
  const auto& projectile=fData->Get(z);
  const auto bracket=FindBracket(projectile,energy);
  if(!std::isfinite(t)) {Fail("Nonfinite Dirac transfer query.");}
  const double top=QueryMaximum(energy,maximum);
  if(t<kTMin||t>top) {return 0.;}
  double result=0.;
  for(int j=0;j<2;++j) {
    const auto& node=*(j==0?bracket.lo:bracket.hi);
    const double weight=j==0?1.-bracket.weight:bracket.weight;
    if(weight==0.) {continue;}
    result+=weight*MappedNode(node,std::min(n,projectile.nmax),top).Density(t);
  }
  return result*ScaleForN(z,n,projectile.nmax);
}

G4double G4CCDiracTable::SampleDelta(G4int z,G4int n,G4double energy,G4double cut,
                                    G4double uNode,G4double uCDF,G4double maximum) const
{
  CheckQuery(n,cut);Validate();
  if(!std::isfinite(uNode)||!std::isfinite(uCDF)||uNode<0.||uNode>1.||uCDF<0.||uCDF>1.) {
    Fail("Dirac sampler requires random numbers in [0,1].");
  }
  const auto& projectile=fData->Get(z);
  const auto bracket=FindBracket(projectile,energy);
  const double top=QueryMaximum(energy,maximum);
  if(cut>=top) {return 0.;}
  const int nn=std::min(n,projectile.nmax);
  const MappedNode low(*bracket.lo,nn,top),high(*bracket.hi,nn,top);
  const double p0=(1.-bracket.weight)*low.Tail(cut)[0];
  const double p1=bracket.weight*high.Tail(cut)[0];
  if(!(p0+p1>0.)) {return 0.;}
  return (p1>0.&&uNode>=p0/(p0+p1)?high:low).Sample(cut,uCDF);
}

G4double G4CCDiracTable::SoftCoefficient(G4int z,G4int zt,G4double energy,G4bool left) const
{
  Validate();return fData->InterpolateStop(z,zt,energy,false,left);
}
G4double G4CCDiracTable::EquilibriumStopping(G4int z,G4int zt,G4double energy,G4bool left) const
{
  Validate();return fData->InterpolateStop(z,zt,energy,true,left);
}
G4double G4CCDiracTable::FullStopping(G4int z,G4int zt,G4int n,G4double energy) const
{
  CheckQuery(n,kT0);
  const double q=static_cast<double>(z-n);
  return SoftCoefficient(z,zt,energy)*q*q+zt*DeltaMoments(z,n,energy,kT0)[1];
}
