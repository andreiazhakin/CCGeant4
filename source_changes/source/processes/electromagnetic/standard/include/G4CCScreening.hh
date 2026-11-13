#ifndef G4CCScreening_hh
#define G4CCScreening_hh 1

#include "globals.hh"

class G4Track;
class G4DynamicParticle;

class G4CCScreening
{
 public:
  static G4double Exponent();
  static void SetExponent(G4double exponent);

  static G4double ScreeningFunction(G4double x);
  static G4double DeltaCharge(G4int Zp, G4int Qe, G4double deltaElectronEnergy);
  static G4double DeltaChargeSquare(G4int Zp, G4int Qe, G4double deltaElectronEnergy);
  static G4double DeltaChargeSquare(const G4DynamicParticle* projectile,
                                    G4double deltaElectronEnergy);

 private:
  static G4double fExponent;
};

#endif
