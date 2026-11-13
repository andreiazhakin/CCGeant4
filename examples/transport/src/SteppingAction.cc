#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include <algorithm>
#include <cmath>

void SteppingAction::UserSteppingAction(const G4Step* step) {
  auto& d = fEvent->Data();
  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  const auto* track = step->GetTrack();
  const bool inTarget = pre->GetPhysicalVolume()->GetLogicalVolume() == fDetector->GetTarget();
  d.edepWorld += step->GetTotalEnergyDeposit();
  if (inTarget) d.edep += step->GetTotalEnergyDeposit();
  if (post->GetStepStatus() == fWorldBoundary) {
    d.escapedEnergy += post->GetKineticEnergy();
    if (track->GetDefinition() == G4Positron::PositronDefinition())
      d.escapedEnergy += 2*electron_mass_c2;
  }
  if (track->GetParentID() == 0 && inTarget) {
    d.finalEnergy = post->GetKineticEnergy();
    d.finalCharge = post->GetCharge()/eplus;
    d.angle = std::acos(std::clamp(d.direction.dot(post->GetMomentumDirection()), -1., 1.));
    d.exited = post->GetStepStatus() == fGeomBoundary;
  }
}
