#pragma once
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"

class RunAction;
struct EventData {
  G4int id = 0;
  G4double initialEnergy = 0, finalEnergy = 0;
  G4double initialCharge = 0, finalCharge = 0;
  G4double edep = 0, edepWorld = 0, escapedEnergy = 0;
  G4double angle = 0;
  G4bool exited = false;
  G4ThreeVector direction{0, 0, 1};
};

class EventAction : public G4UserEventAction {
public:
  explicit EventAction(RunAction* run) : fRun(run) {}
  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;
  EventData& Data() { return fData; }
private:
  RunAction* fRun;
  EventData fData;
};
