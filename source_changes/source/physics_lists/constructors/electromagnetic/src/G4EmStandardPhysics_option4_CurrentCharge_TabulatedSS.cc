#include "G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS.hh"

#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS);

G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS::
G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS(G4int ver, const G4String& name)
  : G4EmStandardPhysics_option4_CurrentCharge(
      ver, IonScatteringMode::TabulatedSingleScattering, name)
{}

namespace {
G4CCPhysicsOptions WithTabulatedSS(G4CCPhysicsOptions options) {
  options.scattering.mode = G4CCScatteringMode::TabulatedSingleScattering;
  return options;
}
}

G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS::
G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS(
  const G4CCPhysicsOptions& options, G4int ver, const G4String& name)
  : G4EmStandardPhysics_option4_CurrentCharge(WithTabulatedSS(options), ver, name)
{}
