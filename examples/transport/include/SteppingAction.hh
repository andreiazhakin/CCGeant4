#pragma once
#include "G4UserSteppingAction.hh"
class DetectorConstruction;
class EventAction;

class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(const DetectorConstruction* detector, EventAction* event)
    : fDetector(detector), fEvent(event) {}
  void UserSteppingAction(const G4Step* step) override;
private:
  const DetectorConstruction* fDetector;
  EventAction* fEvent;
};
