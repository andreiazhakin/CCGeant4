#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(const DetectorConstruction* detector)
  : fDetector(detector) {
  fGun.SetParticleDefinition(G4Geantino::GeantinoDefinition());
  fGun.SetParticleMomentumDirection({0, 0, 1});
}
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  fGun.SetParticlePosition({0, 0, -fDetector->GetThickness()/2 - 1*um});
  fGun.GeneratePrimaryVertex(event);
}
