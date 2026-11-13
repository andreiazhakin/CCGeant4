#ifndef G4CCMottCorrection_h
#define G4CCMottCorrection_h 1

#include "globals.hh"

class G4CCMottCorrection
{
 public:
  static G4bool IsAvailable(G4int zProjectile);

  static G4double RMott(G4int zProjectile, G4double beta,
                        G4double deltaEnergy, G4double tMax);

  static G4bool PolynomialCoefficients(G4int zProjectile, G4double beta,
                                       G4double coeff[5]);

  static G4double FullStoppingCorrection(G4int zProjectile, G4double beta);

  static G4double IntegratedRMottOverT2(G4int zProjectile, G4double beta,
                                        G4double tMin, G4double tMax);

  static G4double G4BornShape(G4double beta2, G4double deltaEnergy,
                              G4double tMax, G4double totalEnergy2,
                              G4double spin);

  static G4double CMottVsG4(G4int zProjectile, G4double beta,
                            G4double deltaEnergy, G4double tMax,
                            G4double beta2, G4double totalEnergy2,
                            G4double spin);

  static G4double MaxRMott(G4int zProjectile, G4double beta,
                           G4double tMin, G4double tMax);

 private:
  static constexpr G4double betaBar = 0.7181228;
};

#endif
