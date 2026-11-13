#include "PhysicsList.hh"

#include "G4CCPhysicsOptions.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option4_CurrentCharge.hh"

#include <cstdlib>
#include <filesystem>
#include <stdexcept>

namespace
{
void SetEnvironment(const char* name, const std::string& value)
{
#ifdef _WIN32
  const int status = _putenv_s(name, value.c_str());
#else
  const int status = setenv(name, value.c_str(), 1);
#endif
  if (status != 0) {
    throw std::runtime_error(std::string("Cannot set ") + name);
  }
}

G4CCScatteringMode ScatteringMode(const std::string& name)
{
  if (name == "msc") return G4CCScatteringMode::Opt4MultipleScattering;
  if (name == "ss") return G4CCScatteringMode::NativeSingleScattering;
  if (name == "tabulated_ss") return G4CCScatteringMode::TabulatedSingleScattering;
  throw std::invalid_argument("Scattering must be msc, ss or tabulated_ss");
}
}

PhysicsList::PhysicsList(const std::string& mode, const std::string& scattering,
                         const std::string& tables, const std::string& elasticTables,
                         double cut)
{
  SetVerboseLevel(0);
  SetDefaultCutValue(cut);

  const auto scatteringMode = ScatteringMode(scattering);
  if (mode == "emopt4") {
    if (scatteringMode != G4CCScatteringMode::Opt4MultipleScattering) {
      throw std::invalid_argument("The emopt4 example uses msc");
    }
    RegisterPhysics(new G4EmStandardPhysics_option4(0));
    return;
  }

  const std::filesystem::path directory(tables);
  const std::string elastic = elasticTables.empty()
    ? (directory / "elastic" / "tables.index.tsv").string() : elasticTables;
  if (scatteringMode == G4CCScatteringMode::TabulatedSingleScattering) {
    const std::filesystem::path path(elastic);
    const auto index = path.filename() == "tables.index.tsv" ? path : path / "tables.index.tsv";
    if (!std::filesystem::is_regular_file(index))
      throw std::runtime_error("Elastic table index not found: " + index.string());
  }

  if (mode == "dirac") {
    G4CCPhysicsOptions cc;
    cc.tables.delta = (directory / "delta").string();
    cc.tables.stopping = (directory / "stopping").string();
    cc.tables.global = (directory / "global.csv").string();
    cc.tables.atimaLow = (directory / "atima_low.csv").string();
    cc.scattering.mode = scatteringMode;
    cc.scattering.tables = elastic;
    RegisterPhysics(new G4EmStandardPhysics_option4_CurrentCharge(cc, 0));
  } else if (mode == "legacy") {
    SetEnvironment("G4CC_GLOBAL_TABLE", (directory / "global.csv").string());
    SetEnvironment("G4CC_ATIMA_Q_DEDX_TABLE", (directory / "legacy_q2.csv").string());
    SetEnvironment("G4CC_ATIMA_EQUILIBRIUM_DEDX_TABLE", (directory / "atima_low.csv").string());
    if (scatteringMode == G4CCScatteringMode::TabulatedSingleScattering) {
      SetEnvironment("G4ELASTIC_TABLE_DIR", elastic);
    }
    RegisterPhysics(new G4EmStandardPhysics_option4_CurrentCharge(0, scatteringMode));
  } else {
    throw std::invalid_argument("Model must be dirac, legacy or emopt4");
  }
}

void PhysicsList::SetCuts()
{
  SetCutsWithDefault();
}
