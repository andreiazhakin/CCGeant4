#ifndef G4CCNonEMChargeResetProcess_hh
#define G4CCNonEMChargeResetProcess_hh 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleChangeForLoss.hh"

class G4CCNonEMChargeResetProcess : public G4VDiscreteProcess
{
 public:
  explicit G4CCNonEMChargeResetProcess(const G4String& name = "CCNonEMChargeReset");
  ~G4CCNonEMChargeResetProcess() override = default;

  G4bool IsApplicable(const G4ParticleDefinition& particle) override;

  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition) override;

  G4double GetMeanFreePath(const G4Track& track,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  G4VParticleChange* PostStepDoIt(const G4Track& track,
                                  const G4Step& step) override;

 private:
  static G4bool ShouldResetAfter(const G4VProcess* process);

  G4ParticleChangeForLoss fParticleChange;
};

#endif
