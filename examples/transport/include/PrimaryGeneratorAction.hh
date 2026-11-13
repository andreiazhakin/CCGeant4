#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  explicit PrimaryGeneratorAction(const DetectorConstruction* detector);
  void GeneratePrimaries(G4Event* event) override;
private:
  const DetectorConstruction* fDetector;
  G4ParticleGun fGun{1};
};
