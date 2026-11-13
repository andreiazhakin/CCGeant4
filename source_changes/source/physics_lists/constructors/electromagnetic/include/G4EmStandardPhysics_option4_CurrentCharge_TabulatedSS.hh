#ifndef G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS_hh
#define G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS_hh 1

#include "G4EmStandardPhysics_option4_CurrentCharge.hh"

class G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS
  : public G4EmStandardPhysics_option4_CurrentCharge
{
 public:
  explicit G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS(
    const G4CCPhysicsOptions& options, G4int ver = 1,
    const G4String& name = "G4EmStandard_opt4_CurrentCharge_Dirac_TabulatedSS");
  explicit G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS(
    G4int ver = 1,
    const G4String& name = "G4EmStandard_opt4_CurrentCharge_TabulatedSS");
  ~G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS() override = default;
};

#endif
