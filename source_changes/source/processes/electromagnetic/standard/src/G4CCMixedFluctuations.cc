#include "G4CCMixedFluctuations.hh"

#include "G4CCCurrentChargeConfig.hh"
#include "G4CCIonChargeState.hh"
#include "G4DynamicParticle.hh"
#include "G4Exception.hh"
#include "G4IonFluctuations.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniversalFluctuation.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>

namespace
{
G4double G4CCAtimaFluctuationCut(const G4MaterialCutsCouple* couple,
                               const G4DynamicParticle* dp,
                               G4double suppliedCut, G4double tmax)
{
  if(couple == nullptr || dp == nullptr || dp->GetDefinition() == nullptr
     || !G4CCCurrentChargeConfig::Instance().UseAtimaEquilibriumDedx()
     || !G4CCIonChargeState::IonisationModelApplies(dp->GetDefinition(),
                                                       couple->GetMaterial())) {
    return suppliedCut;
  }
  const G4double energyPerU = G4CCIonChargeState::KineticEnergyPerNucleon(
    dp->GetDefinition(), dp->GetKineticEnergy())/CLHEP::MeV;
  if(!G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(energyPerU)) {
    return suppliedCut;
  }

  const auto* cuts = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(1);
  const G4int index = couple->GetIndex();
  if(cuts == nullptr || index < 0 || static_cast<std::size_t>(index) >= cuts->size()
     || !std::isfinite((*cuts)[static_cast<std::size_t>(index)])
     || (*cuts)[static_cast<std::size_t>(index)] <= 0.0) {
    G4Exception("G4CCMixedFluctuations::SampleFluctuations", "CCAtimaFluct001",
      FatalException, "Missing or invalid electron production cut for the current couple.");
    return suppliedCut;
  }
  return std::min(tmax, (*cuts)[static_cast<std::size_t>(index)]);
}
}

G4CCMixedFluctuations::G4CCMixedFluctuations(const G4String& name)
  : G4VEmFluctuationModel(name),
    fUniversal(std::make_unique<G4UniversalFluctuation>("CCUniversalFluc")),
    fIon(std::make_unique<G4IonFluctuations>("CCIonFluc"))
{}

G4CCMixedFluctuations::~G4CCMixedFluctuations() = default;

void G4CCMixedFluctuations::InitialiseMe(const G4ParticleDefinition* particle)
{
  fUniversal->InitialiseMe(particle);
  fIon->InitialiseMe(particle);
}

void G4CCMixedFluctuations::SetParticleAndCharge(const G4ParticleDefinition* particle,
                                                 G4double q2)
{
  const G4double q2Safe = std::max(0.0, q2);
  fUniversal->SetParticleAndCharge(particle, q2Safe);
  fIon->SetParticleAndCharge(particle, q2Safe);
}

G4double G4CCMixedFluctuations::IonProbability(const G4DynamicParticle* dp,
                                                   const G4Material* material) const
{
  if(dp == nullptr || dp->GetDefinition() == nullptr) { return 1.0; }

  const auto& mode = G4CCCurrentChargeConfig::Instance().FluctuationMode();
  if(mode == "ion") { return 1.0; }
  if(mode == "universal") { return 0.0; }

  if(material != nullptr &&
     !G4CCIonChargeState::IonisationModelApplies(dp->GetDefinition(), material)) {
    return 1.0;
  }

  const G4double eMeVPerU =
    G4CCIonChargeState::KineticEnergyPerNucleon(dp->GetDefinition(), dp->GetKineticEnergy())/CLHEP::MeV;
  return G4CCCurrentChargeConfig::IonFluctuationMixProbability(eMeVPerU);
}

G4double G4CCMixedFluctuations::SampleFluctuations(
  const G4MaterialCutsCouple* couple,
  const G4DynamicParticle* dp,
  const G4double tcut,
  const G4double tmax,
  const G4double length,
  const G4double meanLoss)
{
  const G4Material* material = (couple == nullptr) ? nullptr : couple->GetMaterial();
  const G4double fluctuationCut = G4CCAtimaFluctuationCut(couple, dp, tcut, tmax);
  const G4double pIon = IonProbability(dp, material);
  if(pIon >= 1.0 || (pIon > 0.0 && G4UniformRand() < pIon)) {
    return fIon->SampleFluctuations(couple, dp, fluctuationCut, tmax, length, meanLoss);
  }
  return fUniversal->SampleFluctuations(couple, dp, fluctuationCut, tmax, length, meanLoss);
}

G4double G4CCMixedFluctuations::Dispersion(
  const G4Material* material,
  const G4DynamicParticle* dp,
  const G4double tcut,
  const G4double tmax,
  const G4double length)
{
  const G4double pIon = IonProbability(dp, material);
  if(pIon <= 0.0) {
    return fUniversal->Dispersion(material, dp, tcut, tmax, length);
  }
  if(pIon >= 1.0) {
    return fIon->Dispersion(material, dp, tcut, tmax, length);
  }
  const G4double dUniversal = fUniversal->Dispersion(material, dp, tcut, tmax, length);
  const G4double dIon = fIon->Dispersion(material, dp, tcut, tmax, length);
  return std::max(0.0, (1.0 - pIon)*dUniversal + pIon*dIon);
}
