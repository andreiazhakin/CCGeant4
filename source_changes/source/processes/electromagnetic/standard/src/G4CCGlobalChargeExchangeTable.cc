#include "G4CCGlobalChargeExchangeTable.hh"

#include "G4CCCurrentChargeConfig.hh"
#include "G4Exception.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <limits>
#include <unordered_set>
#include <sstream>
#include <stdexcept>

namespace
{
std::vector<G4String> SplitCSVLine(const G4String& line)
{
  std::vector<G4String> out;
  G4String cell;
  bool inQuotes = false;
  for(char ch : line) {
    if(ch == '"') {
      inQuotes = !inQuotes;
    } else if(ch == ',' && !inQuotes) {
      out.push_back(cell);
      cell.clear();
    } else {
      cell.push_back(ch);
    }
  }
  out.push_back(cell);
  return out;
}

[[noreturn]] void InvalidGlobal(const G4String& reason)
{
  G4Exception("G4CCGlobalChargeExchangeTable", "CCGlobalData", FatalException, reason.c_str());
  throw std::runtime_error(reason);
}

double StrictNumber(const G4String& text)
{
  std::size_t used=0;
  const double value=std::stod(text,&used);
  if(used!=text.size() || !std::isfinite(value))
    throw std::runtime_error("nonfinite or malformed number");
  return value;
}

int StrictInteger(const G4String& text)
{
  const double value=StrictNumber(text);
  if(std::floor(value)!=value || value<std::numeric_limits<int>::min() ||
     value>std::numeric_limits<int>::max())
    throw std::runtime_error("noninteger or out-of-range integer field");
  return static_cast<int>(value);
}

G4int ToIntRounded(const G4String& s)
{
  return static_cast<G4int>(std::lround(std::stod(s)));
}

G4double ToDoubleOrZero(const G4String& s)
{
  if(s.empty()) { return 0.0; }
  return std::stod(s);
}

G4int RequiredColumn(const std::map<G4String, G4int>& columns,
                     const G4String& name)
{
  const auto it = columns.find(name);
  if(it == columns.end()) {
    G4ExceptionDescription ed;
    ed << "GLOBAL charge-exchange CSV is missing required column '" << name << "'.";
    G4Exception("G4CCGlobalChargeExchangeTable::Load", "G4CC001",
                FatalException, ed);
    throw std::runtime_error("Missing required GLOBAL column: "+name);
  }
  return it->second;
}
}

G4CCGlobalChargeExchangeTable& G4CCGlobalChargeExchangeTable::Instance()
{
  static G4CCGlobalChargeExchangeTable table;
  return table;
}

std::uint32_t G4CCGlobalChargeExchangeTable::PackGroupKey(G4int Zp,
                                                            G4int Zt,
                                                            G4int QeIn)
{
  return (static_cast<std::uint32_t>(Zp & 0xFF) << 16) |
         (static_cast<std::uint32_t>(Zt & 0xFF) << 8) |
          static_cast<std::uint32_t>(QeIn & 0xFF);
}

G4String G4CCGlobalChargeExchangeTable::Trim(const G4String& text)
{
  auto begin = text.begin();
  while(begin != text.end() && std::isspace(static_cast<unsigned char>(*begin))) { ++begin; }
  if(begin == text.end()) { return {}; }
  auto end = text.end();
  do { --end; } while(end != begin && std::isspace(static_cast<unsigned char>(*end)));
  if(std::isspace(static_cast<unsigned char>(*end))) { return {}; }
  return G4String(begin, end + 1);
}

void G4CCGlobalChargeExchangeTable::SetPath(const G4String& path)
{
  std::lock_guard<std::mutex> lock(fLoadMutex);
  fPath = path;
  fLoaded.store(false, std::memory_order_release);
  fLoadAttempted.store(false, std::memory_order_release);
  fData.clear();
  fCompletePairs.clear();
  fStrict=false;
  fNumberOfCurves.store(0, std::memory_order_relaxed);
  fNumberOfPoints.store(0, std::memory_order_relaxed);
}

const G4String& G4CCGlobalChargeExchangeTable::Path() const
{
  return fPath;
}

G4bool G4CCGlobalChargeExchangeTable::IsLoaded() const
{
  return fLoaded.load(std::memory_order_acquire);
}

void G4CCGlobalChargeExchangeTable::EnsureLoaded()
{
  if(fLoaded.load(std::memory_order_acquire) ||
     fLoadAttempted.load(std::memory_order_acquire)) {
    return;
  }

  std::lock_guard<std::mutex> lock(fLoadMutex);
  if(fLoaded.load(std::memory_order_acquire) ||
     fLoadAttempted.load(std::memory_order_acquire)) {
    return;
  }

  LoadUnlocked();
  fLoadAttempted.store(true, std::memory_order_release);
}

void G4CCGlobalChargeExchangeTable::LoadUnlocked()
{
  fStrict=G4CCCurrentChargeConfig::Instance().IsDiracMode();
  fData.clear();fCompletePairs.clear();
  fNumberOfCurves.store(0,std::memory_order_relaxed);
  fNumberOfPoints.store(0,std::memory_order_relaxed);
  if(fPath.empty() && !fStrict) {
    const char* env=std::getenv("G4CC_GLOBAL_TABLE");
    if(env)fPath=env;
  }
  auto unavailable=[&](const char* code,const G4String& text) {
    if(fStrict)InvalidGlobal(text);
    G4Exception("G4CCGlobalChargeExchangeTable::Load",code,JustWarning,text.c_str());
  };
  if(fPath.empty()) {unavailable("G4CC002","No GLOBAL CSV path configured.");return;}
  std::ifstream input(fPath.c_str());
  if(!input) {unavailable("G4CC003","Cannot open GLOBAL CSV: "+fPath);return;}
  G4String header;
  if(!std::getline(input,header)) {unavailable("G4CC004","Empty GLOBAL CSV: "+fPath);return;}
  if(header.size()>=3 && static_cast<unsigned char>(header[0])==0xef &&
     static_cast<unsigned char>(header[1])==0xbb && static_cast<unsigned char>(header[2])==0xbf)
    header.erase(0,3);
  const auto headerCells=SplitCSVLine(header);
  std::map<G4String,G4int> columns;
  for(std::size_t i=0;i<headerCells.size();++i) {
    const auto name=Trim(headerCells[i]);
    if(!columns.emplace(name,static_cast<int>(i)).second && fStrict)
      InvalidGlobal("Duplicate GLOBAL header column: "+name);
  }
  const int cZp=RequiredColumn(columns,"Zp"),cZt=RequiredColumn(columns,"Zt"),
    cE=RequiredColumn(columns,"E_MeV_u"),cQeIn=RequiredColumn(columns,"Qe_in"),
    cQeOut=RequiredColumn(columns,"Qe_out"),cdQe=RequiredColumn(columns,"dQe"),
    cSigma=RequiredColumn(columns,"sigma_cm2");
  const int maxColumn=std::max({cZp,cZt,cE,cQeIn,cQeOut,cdQe,cSigma});
  fData.max_load_factor(0.75F);fData.reserve(500000);
  std::unordered_set<std::uint32_t> pairs;
  G4String line;std::size_t skipped=0,lineNumber=1;
  while(std::getline(input,line)) {
    ++lineNumber;
    if(Trim(line).empty())continue;
    const auto cells=SplitCSVLine(line);
    if(static_cast<int>(cells.size())<=maxColumn ||
       (fStrict && (cells.size()!=headerCells.size() || std::count(line.begin(),line.end(),'"')%2))) {
      if(fStrict)InvalidGlobal("Wrong GLOBAL column/quote count at line "+std::to_string(lineNumber));
      ++skipped;continue;
    }
    try {
      auto integer=[&](int c){return fStrict?StrictInteger(Trim(cells[c])):ToIntRounded(Trim(cells[c]));};
      auto number=[&](int c){return fStrict?StrictNumber(Trim(cells[c])):ToDoubleOrZero(Trim(cells[c]));};
      const int zp=integer(cZp),zt=integer(cZt),n=integer(cQeIn),nout=integer(cQeOut),dn=integer(cdQe);
      const double energy=number(cE);double sigma=number(cSigma);
      if(fStrict) {
        if(zp<29||zp>96||zt<1||zt>96||(zt>=85&&zt<=87)||
           energy<30.||energy>2000.||n<0||n>27||nout<0||nout>27||
           nout!=n+dn||dn==0||std::abs(dn)>2||sigma<0.)
          throw std::runtime_error("row outside supported Z/E/N domain or inconsistent charge transition");
        for(const auto* field:{"Ap","At","GLOBAL_E_MeV_u"}) {
          auto c=columns.find(field);
          if(c!=columns.end() && !(StrictNumber(Trim(cells[c->second]))>0.))
            throw std::runtime_error("invalid reference mass/energy metadata");
        }
        const auto kind=columns.find("kind");
        if(kind!=columns.end() && Trim(cells[kind->second])!=(dn>0?"capture":"loss"))
          throw std::runtime_error("kind disagrees with dQe");
        if(!std::isfinite(static_cast<G4float>(sigma)))
          throw std::runtime_error("sigma exceeds finite runtime storage range");
      } else {
        if(sigma<0.)sigma=0.;
        if(zp<=0||zt<=0||energy<=0.||n<0){++skipped;continue;}
      }
      auto& curves=fData[PackGroupKey(zp,zt,n)];
      auto curve=std::find_if(curves.begin(),curves.end(),[&](const Curve& c){return c.dQe==dn;});
      if(curve==curves.end()) {
        curves.push_back(Curve{dn,nout,{}});curve=curves.end()-1;
        fNumberOfCurves.fetch_add(1,std::memory_order_relaxed);
      } else if(fStrict && curve->QeOut!=nout) {
        throw std::runtime_error("inconsistent Qe_out within a transition curve");
      }
      curve->points.push_back(Point{static_cast<G4float>(energy),nout,static_cast<G4float>(sigma)});
      fNumberOfPoints.fetch_add(1,std::memory_order_relaxed);
      pairs.insert(PackGroupKey(zp,zt,0));
    } catch(const std::exception& error) {
      if(fStrict)InvalidGlobal("Invalid GLOBAL CSV '"+fPath+"', line "+std::to_string(lineNumber)+": "+error.what());
      ++skipped;
    }
  }
  if(input.bad())InvalidGlobal("I/O error while reading GLOBAL CSV: "+fPath);
  for(auto& group:fData)for(auto& curve:group.second) {
    auto& points=curve.points;
    std::sort(points.begin(),points.end(),[](const Point& a,const Point& b){return a.energyMeVPerU<b.energyMeVPerU;});
  }
  if(fStrict) {
    if(fData.empty())InvalidGlobal("GLOBAL CSV contains no valid pair: "+fPath);
    for(const auto pair:pairs) {
      ValidatePairUnlocked(static_cast<int>((pair>>16)&0xff),static_cast<int>((pair>>8)&0xff));
      fCompletePairs.insert(pair);
    }
  }
  const bool loaded=!fData.empty();
  fLoaded.store(loaded,std::memory_order_release);
  if(!loaded)unavailable("G4CC005","No usable GLOBAL rows in "+fPath+"; skipped="+std::to_string(skipped));
}

void G4CCGlobalChargeExchangeTable::ValidatePairUnlocked(G4int zp,G4int zt) const
{
  const G4String label="GLOBAL Zp="+std::to_string(zp)+", Zt="+std::to_string(zt)+": ";
  const std::vector<Point>* reference=nullptr;
  for(int n=0;n<=27;++n) {
    const auto group=fData.find(PackGroupKey(zp,zt,n));
    if(group==fData.end())InvalidGlobal(label+"missing electron state N="+std::to_string(n));
    for(int dn:{-2,-1,1,2}) {
      if(n+dn<0||n+dn>27)continue;
      const auto curve=std::find_if(group->second.begin(),group->second.end(),[&](const Curve& c){return c.dQe==dn;});
      if(curve==group->second.end())InvalidGlobal(label+"missing transition N="+std::to_string(n)+", dN="+std::to_string(dn)+" (zero curves must be explicit)");
      const auto& points=curve->points;
      if(points.size()<2||points.front().energyMeVPerU!=30.F||points.back().energyMeVPerU!=2000.F)
        InvalidGlobal(label+"transition must span 30--2000 MeV/u with at least two nodes");
      for(std::size_t i=1;i<points.size();++i)
        if(!(points[i].energyMeVPerU>points[i-1].energyMeVPerU))
          InvalidGlobal(label+"duplicate energy node in a transition curve");
      if(!reference)reference=&points;
      else {
        if(reference->size()!=points.size())InvalidGlobal(label+"inconsistent energy grids; possible missing row");
        for(std::size_t i=0;i<points.size();++i)
          if((*reference)[i].energyMeVPerU!=points[i].energyMeVPerU)
            InvalidGlobal(label+"inconsistent energy grids; possible missing row");
      }
    }
  }
}

G4bool G4CCGlobalChargeExchangeTable::HasPair(G4int zp,G4int zt) const
{
  return fLoaded.load(std::memory_order_acquire) &&
    fCompletePairs.find(PackGroupKey(zp,zt,0))!=fCompletePairs.end();
}

void G4CCGlobalChargeExchangeTable::RequirePair(G4int zp,G4int zt)
{
  EnsureLoaded();
  if(!HasPair(zp,zt))InvalidGlobal("Missing complete GLOBAL pair Zp="+std::to_string(zp)+", Zt="+std::to_string(zt));
}

void G4CCGlobalChargeExchangeTable::CheckExplicitQuery(G4int zp,G4int zt,G4int n,G4double energy) const
{
  if(!std::isfinite(energy)||energy<=0.||zp<29||zp>92||zt<1||zt>96||
     (zt>=85&&zt<=87)||n<0||n>27)
    InvalidGlobal("Invalid explicit GLOBAL query (Zp=29..92, Zt=1..96 excluding85..87, N=0..27 required)");
}

G4double G4CCGlobalChargeExchangeTable::Interpolate(const std::vector<Point>& points,
                                                    G4double energy)
{
  if(points.empty() || energy <= 0.0) { return 0.0; }

  if(energy < points.front().energyMeVPerU || energy > points.back().energyMeVPerU) {
    return 0.0;
  }

  auto upper = std::upper_bound(points.begin(), points.end(), energy,
                                [](G4double e, const Point& p) {
                                  return e < p.energyMeVPerU;
                                });
  if(upper == points.begin()) { return std::max(0.0, static_cast<G4double>(upper->sigmaCm2)); }
  if(upper == points.end()) { return std::max(0.0, static_cast<G4double>(points.back().sigmaCm2)); }

  const auto lower = upper - 1;
  const G4double e0 = lower->energyMeVPerU;
  const G4double e1 = upper->energyMeVPerU;
  if(e1 <= e0 || e0 <= 0.0) { return std::max(0.0, static_cast<G4double>(lower->sigmaCm2)); }

  const G4double frac = (std::log(energy) - std::log(e0))/(std::log(e1) - std::log(e0));
  const G4double s0 = std::max(0.0, static_cast<G4double>(lower->sigmaCm2));
  const G4double s1 = std::max(0.0, static_cast<G4double>(upper->sigmaCm2));

  if(s0 > 0.0 && s1 > 0.0) {
    return std::exp(std::log(s0) + frac*(std::log(s1) - std::log(s0)));
  }
  return std::max(0.0, s0 + frac*(s1 - s0));
}

std::vector<G4CCGlobalChargeExchangeTable::Transition>
G4CCGlobalChargeExchangeTable::FindTransitions(G4int Zp,
                                               G4int A,
                                               G4int Zt,
                                               G4double kineticEnergyMeVPerU,
                                               G4int QeIn)
{
  std::vector<Transition> out;
  FindTransitions(Zp, A, Zt, kineticEnergyMeVPerU, QeIn, out);
  return out;
}

void G4CCGlobalChargeExchangeTable::FindTransitions(G4int Zp,
                                                    G4int,
                                                    G4int Zt,
                                                    G4double kineticEnergyMeVPerU,
                                                    G4int QeIn,
                                                    std::vector<Transition>& out)
{
  out.clear();
  const bool strict=G4CCCurrentChargeConfig::Instance().IsDiracMode();
  if(strict)CheckExplicitQuery(Zp,Zt,QeIn,kineticEnergyMeVPerU);

  const G4double w = G4CCCurrentChargeConfig::GlobalTransitionWeight(kineticEnergyMeVPerU);
  if(w <= 0.0) { return; }

  const G4double lookupEnergyMeVPerU = kineticEnergyMeVPerU;

  EnsureLoaded();
  if(!fLoaded.load(std::memory_order_acquire)) { return; }

  if(strict)RequirePair(Zp,Zt);
  const auto it = fData.find(PackGroupKey(Zp, Zt, QeIn));
  if(it == fData.end()) { return; }

  out.reserve(it->second.size());
  for(const auto& curve : it->second) {
    const G4double sigma = Interpolate(curve.points, lookupEnergyMeVPerU);
    if(sigma <= 0.0) { continue; }
    const G4int qout = (curve.QeOut >= 0) ? curve.QeOut : (QeIn + curve.dQe);
    out.push_back(Transition{Zt, QeIn, qout, curve.dQe, sigma});
  }
}

G4double G4CCGlobalChargeExchangeTable::InterpolatedSigmaCm2(G4int Zp,
                                                             G4int,
                                                             G4int Zt,
                                                             G4int QeIn,
                                                             G4int dQe,
                                                             G4double kineticEnergyMeVPerU)
{
  const bool strict=G4CCCurrentChargeConfig::Instance().IsDiracMode();
  if(strict)CheckExplicitQuery(Zp,Zt,QeIn,kineticEnergyMeVPerU);
  if(strict && (dQe==0||std::abs(dQe)>2||QeIn+dQe<0||QeIn+dQe>27))
    InvalidGlobal("Invalid explicit GLOBAL transition query");
  const G4double w = G4CCCurrentChargeConfig::GlobalTransitionWeight(kineticEnergyMeVPerU);
  if(w <= 0.0) { return 0.0; }

  const G4double lookupEnergyMeVPerU = kineticEnergyMeVPerU;

  EnsureLoaded();
  if(!fLoaded.load(std::memory_order_acquire)) { return 0.0; }
  if(strict)RequirePair(Zp,Zt);
  const auto it = fData.find(PackGroupKey(Zp, Zt, QeIn));
  if(it == fData.end()) { return 0.0; }
  const auto curve = std::find_if(it->second.begin(), it->second.end(),
                                  [dQe](const Curve& c) { return c.dQe == dQe; });
  if(curve == it->second.end()) { return 0.0; }
  return Interpolate(curve->points, lookupEnergyMeVPerU);
}

std::size_t G4CCGlobalChargeExchangeTable::NumberOfCurves() const
{
  return fNumberOfCurves.load(std::memory_order_acquire);
}

std::size_t G4CCGlobalChargeExchangeTable::NumberOfPoints() const
{
  return fNumberOfPoints.load(std::memory_order_acquire);
}
