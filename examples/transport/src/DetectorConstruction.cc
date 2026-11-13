#include "DetectorConstruction.hh"
#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
#include <stdexcept>

DetectorConstruction::DetectorConstruction() : fThickness(100*um) {
  fMessenger = std::make_unique<G4GenericMessenger>(this, "/example/detector/", "Target");
  fMessenger->DeclareMethod("material", &DetectorConstruction::SetMaterial)
    .SetStates(G4State_PreInit);
  fMessenger->DeclareMethodWithUnit("thickness", "um", &DetectorConstruction::SetThickness)
    .SetStates(G4State_PreInit);
}
DetectorConstruction::~DetectorConstruction() = default;
void DetectorConstruction::SetMaterial(const G4String& value) { fMaterial = value; }
void DetectorConstruction::SetThickness(G4double value) {
  if (!std::isfinite(value) || !(value > 0))
    throw std::invalid_argument("Target thickness must be finite and positive");
  fThickness = value;
}
G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();
  auto* material = nist->FindOrBuildMaterial(fMaterial);
  if (!material) throw std::invalid_argument("Unknown target material: " + fMaterial);
  auto* worldLV = new G4LogicalVolume(
    new G4Box("World", 100*mm, 100*mm, fThickness/2 + 100*mm),
    nist->FindOrBuildMaterial("G4_Galactic"), "World");
  auto* world = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0);
  fTarget = new G4LogicalVolume(new G4Box("Target", 16*mm, 16*mm, fThickness/2),
    material, "Target");
  new G4PVPlacement(nullptr, {}, fTarget, "Target", worldLV, false, 0);
  return world;
}
