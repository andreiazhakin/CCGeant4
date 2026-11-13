#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"

void EventAction::BeginOfEventAction(const G4Event* event) {
  fData = {};
  fData.id = event->GetEventID();
  const auto* particle = event->GetPrimaryVertex()->GetPrimary();
  fData.initialEnergy = particle->GetKineticEnergy();
  fData.initialCharge = particle->GetCharge()/eplus;
  fData.finalEnergy = fData.initialEnergy;
  fData.finalCharge = fData.initialCharge;
  fData.direction = particle->GetMomentumDirection();
}
void EventAction::EndOfEventAction(const G4Event*) { fRun->Write(fData); }
