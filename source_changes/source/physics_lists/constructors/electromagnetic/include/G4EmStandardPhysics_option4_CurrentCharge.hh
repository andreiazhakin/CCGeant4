#ifndef G4EmStandardPhysics_option4_CurrentCharge_hh
#define G4EmStandardPhysics_option4_CurrentCharge_hh 1

#include "G4EmStandardPhysics_option4.hh"
#include "G4CCPhysicsOptions.hh"

class G4EmStandardPhysics_option4_CurrentCharge : public G4EmStandardPhysics_option4
{
 public:
  using IonScatteringMode = G4CCScatteringMode;

  explicit G4EmStandardPhysics_option4_CurrentCharge(
    const G4CCPhysicsOptions& options, G4int ver = 1,
    const G4String& name = "G4EmStandard_opt4_CurrentCharge_Dirac");
  explicit G4EmStandardPhysics_option4_CurrentCharge(
    G4int ver = 1,
    const G4String& name = "G4EmStandard_opt4_CurrentCharge");
  G4EmStandardPhysics_option4_CurrentCharge(
    G4int ver, IonScatteringMode mode,
    const G4String& name = "G4EmStandard_opt4_CurrentCharge");
  ~G4EmStandardPhysics_option4_CurrentCharge() override = default;

  void ConstructProcess() override;
  IonScatteringMode GetIonScatteringMode() const { return ionScatteringMode_; }

 private:
  const IonScatteringMode ionScatteringMode_;
  void ReplaceIonIonisationAndAttachCurrentChargeProcesses();
  void AttachIonSingleScattering();
};

#endif
