#include "RunAction.hh"
#include "EventAction.hh"
#include "G4SystemOfUnits.hh"
#include <filesystem>
#include <iomanip>
#include <stdexcept>

void RunAction::BeginOfRunAction(const G4Run*) {
  if (!IsMaster()) return;
  const std::filesystem::path path(fOutput->filename);
  if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
  fOutput->stream.open(path);
  if (!fOutput->stream) throw std::runtime_error("Cannot open " + fOutput->filename);
  fOutput->stream << std::setprecision(17)
    << "event_id,initial_energy_MeV,final_energy_MeV,loss_MeV,edep_MeV,"
       "initial_charge,final_charge,exit_angle_rad,exited,energy_balance_MeV\n";
}
void RunAction::EndOfRunAction(const G4Run*) {
  if (!IsMaster()) return;
  fOutput->stream.close();
  if (!fOutput->stream) throw std::runtime_error("Cannot write " + fOutput->filename);
}
void RunAction::Write(const EventData& d) {
  std::lock_guard<std::mutex> lock(fOutput->mutex);
  fOutput->stream << d.id << ',' << d.initialEnergy/MeV << ',' << d.finalEnergy/MeV << ','
    << (d.initialEnergy-d.finalEnergy)/MeV << ',' << d.edep/MeV << ','
    << d.initialCharge << ',' << d.finalCharge << ',' << d.angle << ',' << d.exited << ','
    << (d.initialEnergy-d.edepWorld-d.escapedEnergy)/MeV << '\n';
}
