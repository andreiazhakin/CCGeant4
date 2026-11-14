#include "G4CCIonChargeExchangeProcess.hh"

#include "G4CCGlobalChargeExchangeTable.hh"
#include "G4CCCurrentChargeConfig.hh"
#include "G4CCIonChargeState.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ForceCondition.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>

G4CCIonChargeExchangeProcess::G4CCIonChargeExchangeProcess(const G4String& name)
  : G4VDiscreteProcess(name, fElectromagnetic)
{
  SetProcessSubType(91);
}

G4bool G4CCIonChargeExchangeProcess::IsApplicable(const G4ParticleDefinition& particle)
{
  return G4CCIonChargeState::IsIonDefinition(&particle);
}

G4double G4CCIonChargeExchangeProcess::BuildWeights(
  const G4Track& track, std::vector<WeightedTransition>* weights) const
{
  if(weights != nullptr) { weights->clear(); }

  const auto* particle = track.GetParticleDefinition();
  const auto* material = track.GetMaterial();
  if(!G4CCIonChargeState::CurrentChargeModelApplies(particle, material)) { return 0.0; }

  const G4int Zp = G4CCIonChargeState::GetZp(particle);
  const G4int A = G4CCIonChargeState::GetA(particle);
  const G4int Qe = G4CCIonChargeState::GetQe(track);
  const G4double eMeVu = G4CCIonChargeState::KineticEnergyPerNucleon(track);

  if(Zp <= 0 || A <= 0 || Qe < 0 || !std::isfinite(eMeVu) || eMeVu <= 30.0*MeV) { return 0.0; }
  if(G4CCCurrentChargeConfig::GlobalTransitionWeight(eMeVu/MeV) <= 0.0) { return 0.0; }

  const auto* elements = material->GetElementVector();
  const auto* atomDensities = material->GetVecNbOfAtomsPerVolume();
  if(elements == nullptr || atomDensities == nullptr) { return 0.0; }

  G4double sigmaTotal = 0.0;
  std::vector<G4CCGlobalChargeExchangeTable::Transition> transitions;
  transitions.reserve(8);
  const std::size_t nElements = material->GetNumberOfElements();
  for(std::size_t i = 0; i < nElements; ++i) {
    const auto* element = (*elements)[i];
    if(element == nullptr) { continue; }
    const G4int Zt = static_cast<G4int>(std::lround(element->GetZ()));
    if(Zt <= 0 || Zt > G4CCIonChargeState::MaxTargetZ()) { return 0.0; }
    const G4double atomDensity = atomDensities[i];
    if(atomDensity <= 0.0) { continue; }

    G4CCGlobalChargeExchangeTable::Instance().FindTransitions(
      Zp, A, Zt, eMeVu/(MeV), Qe, transitions);
    for(const auto& tr : transitions) {
      const G4double sigmaArea = std::max(0.0, tr.sigmaCm2)*cm2;
      if(sigmaArea <= 0.0) { continue; }
      const G4double weight = atomDensity*sigmaArea;
      sigmaTotal += weight;
      if(weights != nullptr) { weights->push_back(WeightedTransition{Zt, tr.dQe, weight}); }
    }
  }
  return sigmaTotal;
}

G4double G4CCIonChargeExchangeProcess::GetMeanFreePath(const G4Track& track,
                                                       G4double,
                                                       G4ForceCondition* condition)
{
  if(condition != nullptr) { *condition = NotForced; }
  const G4double sigmaTotal = BuildWeights(track, nullptr);
  if(sigmaTotal <= 0.0) { return DBL_MAX; }
  return 1.0/sigmaTotal;
}

G4VParticleChange* G4CCIonChargeExchangeProcess::PostStepDoIt(const G4Track& track,
                                                              const G4Step&)
{
  ClearNumberOfInteractionLengthLeft();

  fParticleChange.InitializeForPostStep(track);

  std::vector<WeightedTransition> weights;
  const G4double total = BuildWeights(track, &weights);
  if(total <= 0.0 || weights.empty()) { return &fParticleChange; }

  G4double sample = G4UniformRand()*total;
  const WeightedTransition* selected = nullptr;
  for(const auto& w : weights) {
    sample -= w.weight;
    if(sample <= 0.0) {
      selected = &w;
      break;
    }
  }
  if(selected == nullptr) { selected = &weights.back(); }

  const G4int Zp = G4CCIonChargeState::GetZp(track.GetParticleDefinition());
  const G4int oldQe = G4CCIonChargeState::GetQe(track);
  const G4int newQe = G4CCIonChargeState::ClampQe(Zp, oldQe + selected->dQe);
  fParticleChange.SetProposedCharge(G4CCIonChargeState::ChargeFromQe(Zp, newQe));

  return &fParticleChange;
}
