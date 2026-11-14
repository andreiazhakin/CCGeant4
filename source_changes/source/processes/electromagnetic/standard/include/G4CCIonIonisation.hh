#ifndef G4CCIonIonisation_hh
#define G4CCIonIonisation_hh 1

#include "G4ionIonisation.hh"

class G4CCIonIonisation : public G4ionIonisation
{
 public:
  explicit G4CCIonIonisation(const G4String& name = "CCIonIoni");
  ~G4CCIonIonisation() override = default;

  void StartTracking(G4Track*) override;

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override;

  G4double PostStepGetPhysicalInteractionLength(
    const G4Track&, G4double previousStepSize, G4ForceCondition*) override;

  G4double AlongStepGetPhysicalInteractionLength(
    const G4Track&, G4double previousStepSize, G4double currentMinimumStep,
    G4double& currentSafety, G4GPILSelection* selection) override;

 protected:
  void InitialiseEnergyLossProcess(const G4ParticleDefinition* particle,
                                   const G4ParticleDefinition* baseParticle) override;

 private:
  G4bool fIsInitialised = false;
  G4bool fPostStepSuppressed = false;
  G4double fBetheBlochThreshold = 0.0;
};

#endif
