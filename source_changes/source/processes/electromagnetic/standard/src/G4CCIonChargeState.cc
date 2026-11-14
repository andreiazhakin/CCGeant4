#include "G4CCIonChargeState.hh"
#include "G4CCDiracPhysics.hh"

#include "G4CCCurrentChargeConfig.hh"

#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4Material.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LossTableManager.hh"
#include "G4EmCorrections.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>

G4bool G4CCIonChargeState::IsIonDefinition(const G4ParticleDefinition* particle)
{
  return particle != nullptr && particle->GetParticleType() == "nucleus" &&
         particle->GetAtomicNumber() > 0 && particle->GetAtomicMass() > 0 &&
         !particle->IsShortLived();
}

G4bool G4CCIonChargeState::IsIonTrack(const G4Track& track)
{
  return IsIonDefinition(track.GetDefinition());
}

G4int G4CCIonChargeState::GetZp(const G4ParticleDefinition* particle)
{
  return IsIonDefinition(particle) ? particle->GetAtomicNumber() : 0;
}

G4int G4CCIonChargeState::GetA(const G4ParticleDefinition* particle)
{
  return IsIonDefinition(particle) ? particle->GetAtomicMass() : 0;
}

G4int G4CCIonChargeState::GetZp(const G4Track& track)
{
  return GetZp(track.GetDefinition());
}

G4int G4CCIonChargeState::GetA(const G4Track& track)
{
  return GetA(track.GetDefinition());
}

G4double G4CCIonChargeState::GetChargeInEplus(const G4DynamicParticle* dynamicParticle)
{
  if(dynamicParticle == nullptr) { return 0.0; }
  return dynamicParticle->GetCharge()/eplus;
}

G4double G4CCIonChargeState::GetChargeInEplus(const G4Track& track)
{
  return GetChargeInEplus(track.GetDynamicParticle());
}

G4int G4CCIonChargeState::GetQe(const G4ParticleDefinition* particle, G4double chargeInEplus)
{
  const G4int Zp = GetZp(particle);
  if(Zp <= 0) { return 0; }
  const auto qRounded = static_cast<G4int>(std::llround(chargeInEplus));
  return ClampQe(Zp, Zp - qRounded);
}

G4int G4CCIonChargeState::GetQe(const G4Track& track)
{
  return GetQe(track.GetDefinition(), GetChargeInEplus(track));
}

G4int G4CCIonChargeState::ClampQe(G4int Zp, G4int Qe)
{
  if(Zp < 0) { Zp = 0; }
  return std::max(0, std::min(Zp, Qe));
}

G4double G4CCIonChargeState::ChargeFromQe(G4int Zp, G4int Qe)
{
  return (static_cast<G4double>(Zp - ClampQe(Zp, Qe))) * eplus;
}

void G4CCIonChargeState::ProposeQe(G4ParticleChangeForLoss& particleChange,
                                    G4int Zp, G4int Qe)
{
  particleChange.SetProposedCharge(ChargeFromQe(Zp, Qe));
}

void G4CCIonChargeState::ProposeBareIon(G4ParticleChangeForLoss& particleChange,
                                         G4int Zp)
{
  particleChange.SetProposedCharge(static_cast<G4double>(std::max(0, Zp))*eplus);
}

G4double G4CCIonChargeState::KineticEnergyPerNucleon(const G4ParticleDefinition* particle,
                                                      G4double kineticEnergy)
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode() && particle && particle->GetPDGMass()>0.)
    return kineticEnergy*(931.494*MeV)/particle->GetPDGMass();
  const G4int A = GetA(particle);
  if(A <= 0) { return kineticEnergy; }
  return kineticEnergy/static_cast<G4double>(A);
}

G4double G4CCIonChargeState::KineticEnergyPerNucleon(const G4Track& track)
{
  return KineticEnergyPerNucleon(track.GetDefinition(), track.GetKineticEnergy());
}

G4int G4CCIonChargeState::MaxTargetZ()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) { return 96; }
  static const G4int maxTargetZ = []() {
    const char* value = std::getenv("G4CC_MAX_TARGET_Z");
    if(value == nullptr || value[0] == '\0') { return 96; }
    const G4int z = std::atoi(value);
    return (z > 0) ? std::min(z, 118) : 96;
  }();
  return maxTargetZ;
}

G4bool G4CCIonChargeState::MaterialIsSupportedTarget(const G4Material* material, G4int maxZt)
{
  if(material == nullptr) { return false; }
  const auto* elements = material->GetElementVector();
  if(elements == nullptr || elements->empty()) { return false; }
  for(const auto* element : *elements) {
    if(element == nullptr) { return false; }
    const G4int z = static_cast<G4int>(std::lround(element->GetZ()));
    if(z <= 0 || z > maxZt) { return false; }
  }
  return true;
}

G4bool G4CCIonChargeState::MaterialIsLightOnly(const G4Material* material, G4int maxZt)
{
  return MaterialIsSupportedTarget(material, maxZt);
}

G4bool G4CCIonChargeState::CurrentChargeModelApplies(const G4Track& track)
{
  return CurrentChargeModelApplies(track.GetDefinition(), track.GetMaterial());
}

G4bool G4CCIonChargeState::CurrentChargeModelApplies(const G4ParticleDefinition* particle,
                                                      const G4Material* material)
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode())
    return G4CCDiracPhysics::Applies(particle,material);
  return IsIonDefinition(particle) && GetZp(particle) >= 29 && MaterialIsSupportedTarget(material, MaxTargetZ());
}

G4double G4CCIonChargeState::CurrentChargeSquare(const G4Track& track)
{
  const G4int Zp = GetZp(track);
  if(Zp <= 0) { return 0.0; }
  const G4int Qe = GetQe(track);
  const G4double q = static_cast<G4double>(Zp - ClampQe(Zp, Qe));
  return q*q;
}

G4double G4CCIonChargeState::EffectiveChargeSquare(const G4Track& track)
{
  const auto* particle = track.GetParticleDefinition();
  const auto* material = track.GetMaterial();
  if(particle == nullptr || material == nullptr) { return CurrentChargeSquare(track); }

  auto* corrections = G4LossTableManager::Instance()->EmCorrections();
  if(corrections == nullptr) { return CurrentChargeSquare(track); }

  corrections->GetParticleCharge(particle, material, track.GetKineticEnergy());
  const G4double q2 = corrections->EffectiveChargeSquareRatio(particle, material,
                                                               track.GetKineticEnergy());
  return (q2 > 0.0) ? q2 : CurrentChargeSquare(track);
}

G4double G4CCIonChargeState::StoppingChargeSquare(const G4Track& track)
{
  if(!CurrentChargeModelApplies(track)) { return EffectiveChargeSquare(track); }

  const G4double eMeVPerU = KineticEnergyPerNucleon(track)/CLHEP::MeV;
  const G4double w = G4CCCurrentChargeConfig::GlobalTransitionWeight(eMeVPerU);
  const G4double q2Current = CurrentChargeSquare(track);
  const G4double q2Effective = EffectiveChargeSquare(track);
  return std::max(0.0, w*q2Current + (1.0 - w)*q2Effective);
}

G4double G4CCIonChargeState::ChargeSquareRatioToBareZ(const G4Track& track)
{
  const G4int Zp = GetZp(track);
  if(Zp <= 0) { return 1.0; }
  const G4double q2 = CurrentChargeSquare(track);
  const G4double z2 = static_cast<G4double>(Zp)*static_cast<G4double>(Zp);
  return q2/z2;
}

G4bool G4CCIonChargeState::IonisationModelApplies(const G4Track& track)
{ return IonisationModelApplies(track.GetDefinition(),track.GetMaterial()); }
G4bool G4CCIonChargeState::IonisationModelApplies(const G4ParticleDefinition* p,const G4Material* m)
{
  return G4CCCurrentChargeConfig::Instance().IsDiracMode()
    ? G4CCDiracPhysics::Applies(p,m) : CurrentChargeModelApplies(p,m);
}
