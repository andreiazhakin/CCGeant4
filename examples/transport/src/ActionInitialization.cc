#include "ActionInitialization.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::BuildForMaster() const {
  SetUserAction(new RunAction(fOutput));
}
void ActionInitialization::Build() const {
  auto* run = new RunAction(fOutput);
  auto* event = new EventAction(run);
  SetUserAction(new PrimaryGeneratorAction(fDetector));
  SetUserAction(run);
  SetUserAction(event);
  SetUserAction(new SteppingAction(fDetector, event));
}
