#include "G4CCIonIonisation.hh"

#include "G4CCLindhardSorensenIonModel.hh"
#include "G4CCCurrentChargeConfig.hh"
#include "G4CCIonChargeState.hh"
#include "G4CCMixedFluctuations.hh"

#include "G4CCBetheBlochModel.hh"
#include "G4CCBraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4BraggIonModel.hh"
#include "G4Electron.hh"
#include "G4Exception.hh"
#include "G4EmParameters.hh"
#include "G4EmProcessSubType.hh"
#include "G4GenericIon.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UniversalFluctuation.hh"
#include "G4IonFluctuations.hh"
#include "G4VEmModel.hh"

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdlib>

namespace
{
G4double G4CCAtimaMaximumFractionalLoss()
{
  static const G4double fraction = []() {
    constexpr G4double defaultFraction = 0.01;
    const char* value = std::getenv("G4CC_ATIMA_MAX_FRACTIONAL_LOSS");
    if(value == nullptr || value[0] == '\0') { return defaultFraction; }
    char* end = nullptr;
    const G4double parsed = std::strtod(value, &end);
    while(end != nullptr && std::isspace(static_cast<unsigned char>(*end))) { ++end; }
    if(end == value || end == nullptr || *end != '\0' || !std::isfinite(parsed)
       || parsed < 1.0e-6 || parsed > 0.02) {
      G4Exception("G4CCIonIonisation", "CCAtimaStep001", FatalException,
        "G4CC_ATIMA_MAX_FRACTIONAL_LOSS must be a finite number in [1e-6, 0.02].");
      return defaultFraction;
    }
    return parsed;
  }();
  return fraction;
}
}

G4CCIonIonisation::G4CCIonIonisation(const G4String& name)
  : G4ionIonisation(name)
{
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(G4Electron::Electron());
}

void G4CCIonIonisation::StartTracking(G4Track* track)
{
  fPostStepSuppressed = false;
  G4ionIonisation::StartTracking(track);
}

G4VParticleChange* G4CCIonIonisation::AlongStepDoIt(
  const G4Track& track, const G4Step& step)
{
  G4CCLindhardSorensenIonModel::TransportScope scope(track);
  return G4ionIonisation::AlongStepDoIt(track, step);
}

G4double G4CCIonIonisation::PostStepGetPhysicalInteractionLength(
  const G4Track& track, G4double previousStepSize, G4ForceCondition* condition)
{
  const G4bool suppress = !HasEnergyLossBiasing() && NumberOfModels() == 1 &&
    dynamic_cast<G4CCLindhardSorensenIonModel*>(GetModelByIndex(0)) != nullptr &&
    G4CCIonChargeState::IonisationModelApplies(track) &&
    G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(
      G4CCIonChargeState::KineticEnergyPerNucleon(track)/MeV);

  if(fPostStepSuppressed && !suppress) {
    ClearNumberOfInteractionLengthLeft();
    currentInteractionLength = DBL_MAX;
    preStepLambda = 0.0;
    mfpKinEnergy = DBL_MAX;
  }

  const G4double baseLimit = G4ionIonisation::PostStepGetPhysicalInteractionLength(
    track, previousStepSize, condition);
  fPostStepSuppressed = suppress;
  if(!suppress) { return baseLimit; }

  currentInteractionLength = DBL_MAX;
  return DBL_MAX;
}

G4double G4CCIonIonisation::AlongStepGetPhysicalInteractionLength(
  const G4Track& track, G4double previousStepSize, G4double currentMinimumStep,
  G4double& currentSafety, G4GPILSelection* selection)
{
  const G4double baseLimit = G4ionIonisation::AlongStepGetPhysicalInteractionLength(
    track, previousStepSize, currentMinimumStep, currentSafety, selection);

  auto* model = dynamic_cast<G4CCLindhardSorensenIonModel*>(EmModel(0));
  if(model == nullptr || !G4CCIonChargeState::IonisationModelApplies(track)
     || !G4CCCurrentChargeConfig::Instance().UseAtimaEquilibriumDedx()) {
    return baseLimit;
  }
  const G4double kineticEnergy = track.GetKineticEnergy();
  const G4double energyPerU = G4CCIonChargeState::KineticEnergyPerNucleon(track)/MeV;
  if(energyPerU < 1.9) { return baseLimit; }

  const G4bool belowAtimaBoundary = energyPerU < 2.0;
  const G4bool boundaryResolutionWindow = energyPerU < 2.0/(1.0 - 0.02);
  const G4double fraction = boundaryResolutionWindow
    ? std::min(G4CCAtimaMaximumFractionalLoss(), 0.00025)
    : G4CCAtimaMaximumFractionalLoss();
  if(!G4CCCurrentChargeConfig::Instance().IsDiracMode() && energyPerU >= 40.0/(1.0 - 0.02)) { return baseLimit; }
  const auto* couple = track.GetMaterialCutsCouple();
  const auto* cuts = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(1);
  const G4int index = couple == nullptr ? -1 : couple->GetIndex();
  if(cuts == nullptr || index < 0 || static_cast<std::size_t>(index) >= cuts->size()) {
    G4Exception("G4CCIonIonisation::AlongStepGetPhysicalInteractionLength",
      "CCAtimaStep002", FatalException, "Missing electron production cut for the current couple.");
    return baseLimit;
  }
  const G4double probeEnergy = !G4CCCurrentChargeConfig::Instance().IsDiracMode() && energyPerU > 40.0
    ? 40.0*MeV*G4CCIonChargeState::GetA(track.GetParticleDefinition()) : kineticEnergy;
  G4double atimaDedx = 0.0;
  if(!belowAtimaBoundary &&
     !model->ComputeAtimaRestrictedDEDX(track.GetMaterial(), track.GetParticleDefinition(),
       probeEnergy, (*cuts)[static_cast<std::size_t>(index)],
       G4CCIonChargeState::GetChargeInEplus(track), atimaDedx)) {
    return baseLimit;
  }

  const G4double seedDedx = GetDEDX(kineticEnergy, couple);
  const G4double referenceRange = GetRange(kineticEnergy, couple);
  if(!std::isfinite(atimaDedx) || atimaDedx < 0.0
     || !std::isfinite(seedDedx) || seedDedx <= 0.0
     || !std::isfinite(referenceRange) || referenceRange <= 0.0) {
    G4ExceptionDescription description;
    description << "Invalid ATIMA/reference step estimate at E/u=" << energyPerU
                << " MeV: ATIMA DEDX=" << atimaDedx << ", seed DEDX=" << seedDedx
                << ", reference range=" << referenceRange;
    G4Exception("G4CCIonIonisation::AlongStepGetPhysicalInteractionLength",
      "CCAtimaStep003", FatalException, description);
    return baseLimit;
  }

  const G4double lossLimit = fraction*kineticEnergy/std::max(atimaDedx, seedDedx);
  return std::min(baseLimit, std::min(lossLimit, 0.5*referenceRange));
}

void G4CCIonIonisation::InitialiseEnergyLossProcess(const G4ParticleDefinition* particle,
                                                    const G4ParticleDefinition* baseParticle)
{
  if(fIsInitialised) { return; }

  const auto* genericIon = G4GenericIon::GenericIon();
  const G4int pdg = particle == nullptr ? 0 : particle->GetPDGEncoding();

  const G4ParticleDefinition* theBaseParticle = nullptr;
  if(particle == baseParticle) {
    theBaseParticle = nullptr;
  } else if(baseParticle != nullptr) {
    theBaseParticle = baseParticle;
  } else if(particle == genericIon || pdg == 1000020040) {
    theBaseParticle = nullptr;
  } else {
    theBaseParticle = genericIon;
  }
  SetBaseParticle(theBaseParticle);

  if(particle != nullptr) {
    fBetheBlochThreshold = 2.0*MeV*particle->GetPDGMass()/proton_mass_c2;
  } else {
    fBetheBlochThreshold = 2.0*MeV;
  }

  auto* param = G4EmParameters::Instance();
  const G4double emin = param->MinKinEnergy();
  G4double emax = param->MaxKinEnergy();

  if(FluctModel() == nullptr) {
    const G4String& mode = G4CCCurrentChargeConfig::Instance().FluctuationMode();
    if(G4CCCurrentChargeConfig::Instance().UseAtimaEquilibriumDedx()) {
      SetFluctModel(new G4CCMixedFluctuations());
    } else if(mode == "ion") {
      SetFluctModel(new G4IonFluctuations());
    } else if(mode == "mix") {
      SetFluctModel(new G4CCMixedFluctuations());
    } else {
      SetFluctModel(new G4UniversalFluctuation());
    }
  }

  if(EmModel(0) == nullptr) {
    if(pdg == 1000020040) {
      SetEmModel(new G4BraggIonModel());
    } else {
      SetEmModel(new G4CCLindhardSorensenIonModel());
    }
  }

  EmModel(0)->SetLowEnergyLimit(emin);
  EmModel(0)->SetHighEnergyLimit(emax);
  AddEmModel(1, EmModel(0), FluctModel());

  fIsInitialised = true;
}
