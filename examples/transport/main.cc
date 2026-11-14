#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"

#include "G4CCPhysicsOptions.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option4_CurrentCharge.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
  try {
    std::string model = "dirac";
    std::string scattering = "msc";
    std::string tables = "data";
    std::string macro = "macros/smoke.mac";
    std::string output = "events.csv";

    for (int i = 1; i < argc; ++i) {
      const std::string key = argv[i];
      if (key == "--help") {
        std::cout << "cc_transport --macro FILE --model dirac|legacy|emopt4\n"
                     "  --scattering msc|ss|tabulated_ss --tables DIRECTORY --output CSV\n";
        return 0;
      }
      if (++i == argc) throw std::invalid_argument("Missing value for " + key);
      const std::string value = argv[i];
      if (key == "--macro") macro = value;
      else if (key == "--model") model = value;
      else if (key == "--scattering") scattering = value;
      else if (key == "--tables") tables = value;
      else if (key == "--output") output = value;
      else throw std::invalid_argument("Unknown option: " + key);
    }

    auto runManager = std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager());
    auto* detector = new DetectorConstruction;
    runManager->SetUserInitialization(detector);

    auto* physics = new G4VModularPhysicsList;
    physics->SetVerboseLevel(0);

    if (model == "emopt4") {
      if (scattering != "msc")
        throw std::invalid_argument("The emopt4 example uses msc");
      physics->RegisterPhysics(new G4EmStandardPhysics_option4(0));
    } else {
      G4CCPhysicsOptions cc;
      cc.tables.global = tables + "/global.csv";
      cc.tables.atimaLow = tables + "/atima_low.csv";

      if (model == "dirac") {
        cc.mode = G4CCMode::Dirac;
        cc.tables.delta = tables + "/delta";
        cc.tables.stopping = tables + "/stopping";
      } else if (model == "legacy") {
        cc.mode = G4CCMode::Legacy;
        cc.tables.legacyStopping = tables + "/legacy_q2.csv";
      } else {
        throw std::invalid_argument("Model must be dirac, legacy or emopt4");
      }

      if (scattering == "msc") {
        cc.scattering.mode = G4CCScatteringMode::Opt4MultipleScattering;
      } else if (scattering == "ss") {
        cc.scattering.mode = G4CCScatteringMode::NativeSingleScattering;
      } else if (scattering == "tabulated_ss") {
        cc.scattering.mode = G4CCScatteringMode::TabulatedSingleScattering;
        cc.scattering.tables = tables + "/elastic/tables.index.tsv";
      } else {
        throw std::invalid_argument("Scattering must be msc, ss or tabulated_ss");
      }

      physics->RegisterPhysics(new G4EmStandardPhysics_option4_CurrentCharge(cc, 0));
    }

    runManager->SetUserInitialization(physics);
    auto runOutput = std::make_shared<RunOutput>(output);
    runManager->SetUserInitialization(new ActionInitialization(detector, runOutput));

    const auto status = G4UImanager::GetUIpointer()->ApplyCommand("/control/execute " + macro);
    return status == 0 ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
