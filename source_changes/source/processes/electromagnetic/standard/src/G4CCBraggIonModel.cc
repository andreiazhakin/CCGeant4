#include "G4CCBraggIonModel.hh"

#include "G4CCIonChargeState.hh"
#include "G4CCCurrentChargeConfig.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VEmFluctuationModel.hh"

#include <algorithm>
#include <cmath>

namespace
{
G4bool G4CCUseLowEnergyVanillaFallback(const G4ParticleDefinition* particle,
                                       const G4Material* material,
                                       G4double kineticEnergy)
{
  if(!G4CCIonChargeState::IonisationModelApplies(particle, material)) { return false; }
  const G4double eMeVPerU =
    G4CCIonChargeState::KineticEnergyPerNucleon(particle, kineticEnergy)/CLHEP::MeV;
  return eMeVPerU <= (G4CCCurrentChargeConfig::Instance().IsDiracMode()
    ? G4CCCurrentChargeConfig::Low0MeVPerU() : 30.0);
}
}

G4CCBraggIonModel::G4CCBraggIonModel(const G4ParticleDefinition* p,
                                     const G4String& name)
  : G4BraggModel(p, name)
{}

G4double G4CCBraggIonModel::ChargeSquareRatio(const G4Track& track)
{
  if(!G4CCIonChargeState::IonisationModelApplies(track) ||
     G4CCUseLowEnergyVanillaFallback(track.GetParticleDefinition(),
                                     track.GetMaterial(),
                                     track.GetKineticEnergy())) {
    fHasCachedCurrentCharge = false;
    return G4BraggModel::ChargeSquareRatio(track);
  }

  fCachedChargeInEplus = G4CCIonChargeState::GetChargeInEplus(track);
  fHasCachedCurrentCharge = true;

  const G4double q2 = std::max(0.0, G4CCIonChargeState::StoppingChargeSquare(track));
  if(auto* fluct = GetModelOfFluctuations()) {
    fluct->SetParticleAndCharge(track.GetParticleDefinition(), q2);
  }
  return q2;
}

G4double G4CCBraggIonModel::CrossSectionPerVolume(const G4Material* material,
                                                  const G4ParticleDefinition* particle,
                                                  G4double kineticEnergy,
                                                  G4double cutEnergy,
                                                  G4double maxEnergy)
{
  if(G4CCIonChargeState::IonisationModelApplies(particle, material) &&
     G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(
       G4CCIonChargeState::KineticEnergyPerNucleon(particle, kineticEnergy)/CLHEP::MeV)) {
    return 0.0;
  }
  return G4BraggModel::CrossSectionPerVolume(material, particle, kineticEnergy,
                                             cutEnergy, maxEnergy);
}

G4double G4CCBraggIonModel::GetChargeSquareRatio(const G4ParticleDefinition* particle,
                                                 const G4Material* material,
                                                 G4double kineticEnergy)
{
  if(G4CCIonChargeState::IonisationModelApplies(particle, material) &&
     !G4CCUseLowEnergyVanillaFallback(particle, material, kineticEnergy)) {
    return 1.0;
  }
  return G4BraggModel::GetChargeSquareRatio(particle, material, kineticEnergy);
}

G4double G4CCBraggIonModel::GetParticleCharge(const G4ParticleDefinition* particle,
                                              const G4Material* material,
                                              G4double kineticEnergy)
{
  if(G4CCIonChargeState::IonisationModelApplies(particle, material) &&
     !G4CCUseLowEnergyVanillaFallback(particle, material, kineticEnergy)) {
    if(fHasCachedCurrentCharge) { return fCachedChargeInEplus*eplus; }
    return static_cast<G4double>(G4CCIonChargeState::GetZp(particle))*eplus;
  }
  return G4BraggModel::GetParticleCharge(particle, material, kineticEnergy);
}

void G4CCBraggIonModel::CorrectionsAlongStep(const G4Material* material,
                                             const G4ParticleDefinition* particle,
                                             const G4double kineticEnergy,
                                             const G4double cutEnergy,
                                             const G4double& length,
                                             G4double& eloss)
{
  if(G4CCIonChargeState::IonisationModelApplies(particle, material) &&
     !G4CCUseLowEnergyVanillaFallback(particle, material, kineticEnergy)) {
    return;
  }
  G4BraggModel::CorrectionsAlongStep(material, particle, kineticEnergy,
                                     cutEnergy, length, eloss);
}

void G4CCBraggIonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* secondaries,
                                          const G4MaterialCutsCouple* couple,
                                          const G4DynamicParticle* dynamicParticle,
                                          G4double tmin,
                                          G4double maxEnergy)
{
  const auto* material = couple == nullptr ? nullptr : couple->GetMaterial();
  if(dynamicParticle != nullptr &&
     G4CCIonChargeState::IonisationModelApplies(dynamicParticle->GetDefinition(), material) &&
     G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(
       G4CCIonChargeState::KineticEnergyPerNucleon(dynamicParticle->GetDefinition(),
                                                 dynamicParticle->GetKineticEnergy())/CLHEP::MeV)) {
    return;
  }
  G4BraggModel::SampleSecondaries(secondaries, couple, dynamicParticle, tmin, maxEnergy);
}
