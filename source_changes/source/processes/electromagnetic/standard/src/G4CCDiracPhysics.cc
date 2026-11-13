#include "G4CCDiracPhysics.hh"
#include "G4CCDiracTable.hh"
#include "G4CCCurrentChargeConfig.hh"
#include "G4CCGlobalChargeExchangeTable.hh"
#include "G4CCIonChargeState.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Exception.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <stdexcept>

namespace G4CCDiracPhysics {
G4CCDiracTable& Tables() {
  auto& table = G4CCDiracTable::Instance();
  static std::once_flag loaded;
  std::call_once(loaded, [&] {
    const auto& cfg = G4CCCurrentChargeConfig::Instance().Options();
    table.Configure(cfg.tables.delta, cfg.tables.stopping);
    table.Validate();
  });
  return table;
}
namespace {
[[noreturn]] void DomainError(const G4String& message) {
  G4Exception("G4CCDiracPhysics","CCDiracDomain",FatalException,message.c_str());
  throw std::runtime_error(message);
}
}
G4bool MaterialInDomain(const G4Material* m) {
  if(!m || !m->GetElementVector() || m->GetElementVector()->empty())return false;
  bool supported=true;
  for(const auto* el:*m->GetElementVector()) {
    if(!el)return false;
    const int zt=el->GetZasInt();
    if(zt>=85&&zt<=87)
      DomainError("Explicit Current Charge has no GLOBAL target data for Zt="+std::to_string(zt));
    if(zt<1||zt>96)supported=false;
  }
  return supported;
}
void RequireMaterialCoverage(G4int z,const G4Material* m) {
  if(z<29||z>92)DomainError("Explicit Current Charge requires projectile Zp=29..92; requested Zp="+std::to_string(z));
  auto& data=Tables();
  if(!data.HasProjectile(z))DomainError("Missing required Dirac projectile data for Zp="+std::to_string(z));
  for(const auto* el:*m->GetElementVector()) {
    const int zt=el->GetZasInt();
    if(!data.HasPair(z,zt))DomainError("Missing required Dirac stopping pair Zp="+std::to_string(z)+", Zt="+std::to_string(zt));
    G4CCGlobalChargeExchangeTable::Instance().RequirePair(z,zt);
  }
}
G4bool Applies(const G4ParticleDefinition* p,const G4Material* m) {
  if(!p||!m||!G4CCIonChargeState::IsIonDefinition(p))return false;
  const int z=G4CCIonChargeState::GetZp(p);
  if(z<29)return false;
  if(!MaterialInDomain(m))return false;
  if(z<=92 && !Tables().HasProjectile(z))return false;
  RequireMaterialCoverage(z,m);
  return true;
}
void CheckCut(G4double cut) {
  if(!std::isfinite(cut) || cut < G4CCDiracTable::MinTransfer()*MeV) {
    G4Exception("G4CCDiracPhysics", "CCDiracCut", FatalException,
      "Dirac delta production requires an electron energy cut >= 0.99 keV in each applicable material/region.");
    throw std::runtime_error("Dirac electron cut below 0.99 keV or nonfinite");
  }
}
namespace {
struct Mapping { double eu, cut, top; };
Mapping Map(const G4ParticleDefinition* p, double kinetic, double cut) {
  CheckCut(cut);
  const double mass=p->GetPDGMass()/MeV, tau=(kinetic/MeV)/mass;
  const double eu=G4CCIonChargeState::KineticEnergyPerNucleon(p,kinetic)/MeV;
  const double r=(electron_mass_c2/MeV)/mass;
  const double top=2*(electron_mass_c2/MeV)*tau*(tau+2)/(1+2*(tau+1)*r+r*r);
  return {eu,cut/MeV,top};
}
}
std::array<G4double,3> Moments(const G4ParticleDefinition* p,int n,double k,double cut) {
  const auto m=Map(p,k,cut);
  if(m.cut>=m.top) return {0.,0.,0.};
  return Tables().DeltaMoments(G4CCIonChargeState::GetZp(p),n,m.eu,m.cut,m.top);
}
G4double Sample(const G4ParticleDefinition* p,int n,double k,double cut,double un,double uc) {
  const auto m=Map(p,k,cut);
  if(m.cut>=m.top) return 0.;
  return Tables().SampleDelta(G4CCIonChargeState::GetZp(p),n,m.eu,m.cut,un,uc,m.top)*MeV;
}
G4double FullStopping(const G4Material* mat,const G4ParticleDefinition* p,int n,double k) {
  const double eu=G4CCIonChargeState::KineticEnergyPerNucleon(p,k)/MeV;
  const auto* elements=mat->GetElementVector();
  const auto* density=mat->GetVecNbOfAtomsPerVolume();
  double value=0.;
  for(std::size_t i=0;i<elements->size();++i)
    value += density[i]*Tables().FullStopping(G4CCIonChargeState::GetZp(p),
               (*elements)[i]->GetZasInt(),n,eu)*MeV*cm2;
  return value;
}
G4double EquilibriumStopping(const G4Material* mat,const G4ParticleDefinition* p,double k) {
  const double eu=G4CCIonChargeState::KineticEnergyPerNucleon(p,k)/MeV;
  const auto* elements=mat->GetElementVector();
  const auto* density=mat->GetVecNbOfAtomsPerVolume();
  double value=0.;
  for(std::size_t i=0;i<elements->size();++i)
    value += density[i]*Tables().EquilibriumStopping(G4CCIonChargeState::GetZp(p),
               (*elements)[i]->GetZasInt(),eu)*MeV*cm2;
  return value;
}
G4double HardStopping(const G4Material* mat,const G4ParticleDefinition* p,int n,double k,double cut) {
  return mat->GetElectronDensity()*Moments(p,n,k,cut)[1]*MeV*cm2;
}
}
