#pragma once
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <memory>

class G4GenericMessenger;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction() override;
  G4VPhysicalVolume* Construct() override;
  G4double GetThickness() const { return fThickness; }
  const G4LogicalVolume* GetTarget() const { return fTarget; }
private:
  void SetMaterial(const G4String& value);
  void SetThickness(G4double value);
  G4String fMaterial = "G4_Si";
  G4double fThickness;
  G4LogicalVolume* fTarget = nullptr;
  std::unique_ptr<G4GenericMessenger> fMessenger;
};
