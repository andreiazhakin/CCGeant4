#ifndef G4CCPhysicsOptions_hh
#define G4CCPhysicsOptions_hh 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

enum class G4CCScatteringMode {
  Opt4MultipleScattering,
  NativeSingleScattering,
  TabulatedSingleScattering
};

struct G4CCPhysicsOptions {
  struct Tables {
    G4String delta;
    G4String stopping;
    G4String atimaLow;
    G4String global;
  } tables;
  struct Transition {
    G4double low0 = 30*CLHEP::MeV;
    G4double low1 = 40*CLHEP::MeV;
    G4double high0 = 1500*CLHEP::MeV;
    G4double high1 = 2000*CLHEP::MeV;
  } transition;
  struct Scattering {
    G4CCScatteringMode mode = G4CCScatteringMode::Opt4MultipleScattering;
    G4String tables;
  } scattering;
  G4String fluctuationMode = "mix";
  G4double fluctuationMixPower = 1.;
};

#endif
