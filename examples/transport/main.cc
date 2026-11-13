#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    std::string model="dirac", scattering="msc", tables="data", elastic;
    std::string macro="macros/smoke.mac", output="events.csv";
    int threads=1;
    long seed=6110101;
    double cut=25*um;
    for (int i=1; i<argc; ++i) {
      const std::string key=argv[i];
      if (key=="--help") {
        std::cout << "cc_transport --model dirac|legacy|emopt4 --scattering msc|ss|tabulated_ss\n"
                     "  --tables DIRECTORY --elastic-tables INDEX --macro FILE --output CSV\n"
                     "  --cut-um 25 --threads 1 --seed 6110101\n";
        return 0;
      }
      if (++i==argc) throw std::invalid_argument("Missing value for " + key);
      const std::string value=argv[i];
      std::size_t parsed=value.size();
      if (key=="--model") model=value;
      else if (key=="--scattering") scattering=value;
      else if (key=="--tables") tables=value;
      else if (key=="--elastic-tables") elastic=value;
      else if (key=="--macro") macro=value;
      else if (key=="--output") output=value;
      else if (key=="--cut-um") cut=std::stod(value, &parsed)*um;
      else if (key=="--threads") threads=std::stoi(value, &parsed);
      else if (key=="--seed") seed=std::stol(value, &parsed);
      else throw std::invalid_argument("Unknown option: " + key);
      if (parsed!=value.size()) throw std::invalid_argument("Invalid value for " + key + ": " + value);
    }
    if (threads<1 || !std::isfinite(cut) || !(cut>0) || seed<=0)
      throw std::invalid_argument("Threads and seed must be positive; cut must be finite and positive");
    G4Random::setTheSeed(seed);
    auto manager=std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager());
    manager->SetNumberOfThreads(threads);
    auto* detector=new DetectorConstruction;
    manager->SetUserInitialization(detector);
    manager->SetUserInitialization(new PhysicsList(model, scattering, tables, elastic, cut));
    manager->SetUserInitialization(new ActionInitialization(detector, std::make_shared<RunOutput>(output)));
    return G4UImanager::GetUIpointer()->ApplyCommand("/control/execute " + macro)==0 ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
