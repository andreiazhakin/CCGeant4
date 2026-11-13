#pragma once
#include "G4VUserActionInitialization.hh"
#include <memory>
#include <utility>
class DetectorConstruction;
struct RunOutput;

class ActionInitialization : public G4VUserActionInitialization {
public:
  ActionInitialization(const DetectorConstruction* detector, std::shared_ptr<RunOutput> output)
    : fDetector(detector), fOutput(std::move(output)) {}
  void BuildForMaster() const override;
  void Build() const override;
private:
  const DetectorConstruction* fDetector;
  std::shared_ptr<RunOutput> fOutput;
};
