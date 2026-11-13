#pragma once
#include "G4UserRunAction.hh"
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <utility>

struct EventData;
struct RunOutput {
  explicit RunOutput(std::string filename) : filename(std::move(filename)) {}
  std::string filename;
  std::ofstream stream;
  std::mutex mutex;
};
class RunAction : public G4UserRunAction {
public:
  explicit RunAction(std::shared_ptr<RunOutput> output) : fOutput(std::move(output)) {}
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;
  void Write(const EventData& data);
private:
  std::shared_ptr<RunOutput> fOutput;
};
