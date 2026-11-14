#ifndef G4CCIonChargeExchangeProcess_hh
#define G4CCIonChargeExchangeProcess_hh 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleChangeForLoss.hh"

#include <vector>

class G4CCIonChargeExchangeProcess : public G4VDiscreteProcess
{
 public:
  explicit G4CCIonChargeExchangeProcess(const G4String& name = "CCIonChargeExchange");
  ~G4CCIonChargeExchangeProcess() override = default;

  G4bool IsApplicable(const G4ParticleDefinition& particle) override;

  G4double GetMeanFreePath(const G4Track& track,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  G4VParticleChange* PostStepDoIt(const G4Track& track,
                                  const G4Step& step) override;

 private:
  struct WeightedTransition
  {
    G4int Zt = 0;
    G4int dQe = 0;
    G4double weight = 0.0;
  };

  G4double BuildWeights(const G4Track& track,
                        std::vector<WeightedTransition>* weights) const;

  G4ParticleChangeForLoss fParticleChange;
};

#endif
