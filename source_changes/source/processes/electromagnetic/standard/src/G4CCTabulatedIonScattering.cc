#include "G4CCTabulatedIonScattering.hh"
#include "G4CCTabulatedIonScatteringModel.hh"
#include "G4IonCoulombScatteringModel.hh"
#include "G4Exception.hh"
#include "G4EmParameters.hh"
#include "G4EmProcessSubType.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Track.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utility>

G4CCTabulatedIonScattering::G4CCTabulatedIonScattering(
    std::shared_ptr<const G4CCElasticTableData> data,
    G4bool useTables, G4double recoilThreshold, G4bool allowEnvironment)
    : G4VEmProcess(useTables ? "CCTabulatedIonScattering" : "CoulombScat"),
      useTables_(useTables) {
  if (!(recoilThreshold >= 0.) || !std::isfinite(recoilThreshold)) {
    G4Exception("G4CCTabulatedIonScattering", "CCElasticConfig001",
                FatalException, "Invalid recoil threshold");
    throw std::invalid_argument("Invalid recoil threshold");
  }
  if (useTables_) {
    try {
      if (!data && allowEnvironment) data = G4CCElasticTableData::LoadFromEnvironment();
      if(!data) {
        static std::once_flag warning;
        std::call_once(warning, []() {
          G4Exception("G4CCTabulatedIonScattering", "CCElasticMissingTables",
                      JustWarning, "No single-scattering table/index was supplied or found; using native Geant4 single scattering.");
        });
      }
      auto* model = new G4CCTabulatedIonScatteringModel(std::move(data));
      model->SetRecoilThreshold(recoilThreshold);
      SetEmModel(model);
    } catch (const std::exception& error) {
      G4Exception("G4CCTabulatedIonScattering", "CCElasticData001",
                  FatalException, error.what());
      throw;
    }
  } else {
    auto* model = new G4IonCoulombScatteringModel();
    model->SetRecoilThreshold(recoilThreshold);
    SetEmModel(model);
  }
  SetBuildTableFlag(false);
  SetStartFromNullFlag(false);
  SetSplineFlag(false);
  SetSecondaryParticle(G4Proton::Proton());
  SetProcessSubType(fCoulombScattering);
}

G4bool G4CCTabulatedIonScattering::IsApplicable(
    const G4ParticleDefinition& particle) {
  return particle.GetPDGCharge() != 0. &&
         (particle.GetParticleType() == "nucleus" ||
          &particle == G4Proton::Proton());
}

void G4CCTabulatedIonScattering::PreparePhysicsTable(
    const G4ParticleDefinition& particle) {
  G4VEmProcess::PreparePhysicsTable(particle);
  if (EmModel()) EmModel()->SetPolarAngleLimit(0.);
}

void G4CCTabulatedIonScattering::StartTracking(G4Track* track) {
  G4VEmProcess::StartTracking(track);
  coupleIdxLambda = std::numeric_limits<std::size_t>::max();
}

G4double G4CCTabulatedIonScattering::MinPrimaryEnergy(
    const G4ParticleDefinition*, const G4Material*) {
  return 0.;
}

void G4CCTabulatedIonScattering::InitialiseProcess(
    const G4ParticleDefinition* particle) {
  if (!isInitialised_) {
    auto* parameters = G4EmParameters::Instance();
    auto* model = EmModel(0);
    model->SetLowEnergyLimit(std::max(parameters->MinKinEnergy(),
                                     model->LowEnergyLimit()));
    model->SetHighEnergyLimit(std::min(parameters->MaxKinEnergy(),
                                      model->HighEnergyLimit()));
    AddEmModel(1, model);
    if (particle->GetParticleName() != "GenericIon") SetVerboseLevel(0);
    isInitialised_ = true;
  }
  SetCrossSectionType(fEmNoIntegral);
  SetBuildTableFlag(false);
  SetStartFromNullFlag(false);
  if (EmModel()) EmModel()->SetPolarAngleLimit(0.);
}

void G4CCTabulatedIonScattering::StreamProcessInfo(std::ostream& out) const {
  out << "      Built-in full ion single scattering: 0 <= theta_CM <= pi; "
         "direct per-step sigma; no integral envelope or MSC split. "
      << (useTables_ ? "Z-pair tables with native fallback."
                     : "Standard G4IonCoulombScatteringModel.")
      << '\n';
}

void G4CCTabulatedIonScattering::ProcessDescription(std::ostream& out) const {
  out << "Built-in full-range ion elastic single scattering. "
         "Table mode replaces the differential measure by a table selected "
         "only by nuclear (Z_projectile,Z_target), evaluated at actual T/A. "
         "Wholly uncovered materials/energies delegate to the standard "
         "G4IonCoulombScatteringModel. Uncovered targets in mixed coverage "
         "use its native Coulomb helper after one target selection. "
         "Target selection, recoil and primary kinetic changes belong to "
         "this one process; charge and isotope identity are preserved.\n";
  G4VEmProcess::ProcessDescription(out);
}
