// The fallback cross-section and angular law are supplied by Geant4's
// G4IonCoulombCrossSection. Target selection and recoil handling below adapt
// G4IonCoulombScatteringModel in this separately named built-in model.
// Geant4: Copyright (c) the Geant4 Collaboration, Geant4 Software License
// https://geant4.web.cern.ch/download/license

#include "G4CCTabulatedIonScatteringModel.hh"
#include "G4CCElasticTableData.hh"

#include "G4Exception.hh"
#include "G4IonCoulombCrossSection.hh"
#include "G4IonCoulombScatteringModel.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4NistManager.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace tabulated {
namespace {
void WarnMissingPair(int zp, int zt) {
  static std::mutex mutex;
  static std::set<std::pair<int,int>> warned;
  bool first = false;
  {
    std::lock_guard<std::mutex> guard(mutex);
    first = warned.emplace(zp,zt).second;
  }
  if(first) {
    const std::string message = "No tabulated single-scattering pair Zp=" +
      std::to_string(zp) + ", Zt=" + std::to_string(zt) +
      "; using native Geant4 single scattering for this pair.";
    G4Exception("G4CCTabulatedIonScatteringModel", "CCElasticMissingPair", JustWarning, message.c_str());
  }
}

[[noreturn]] void Fail(const char* where, const std::string& message) {
  G4Exception(where, "TabulatedElastic001", FatalException, message.c_str());
  throw std::runtime_error(message);
}

struct NucleusKey {
  const G4ParticleDefinition* projectile;
  int z, a, isomer;
  bool operator==(const NucleusKey& other) const {
    return projectile == other.projectile && z == other.z && a == other.a &&
           isomer == other.isomer;
  }
};
struct NucleusHash {
  std::size_t operator()(const NucleusKey& k) const {
    auto h = std::hash<const void*>{}(k.projectile);
    h ^= static_cast<std::size_t>(k.z * 1000003 + k.a * 1009 + k.isomer);
    return h;
  }
};
struct MaterialKey {
  const G4Material* material;
  const G4ParticleDefinition* projectile;
  bool operator==(const MaterialKey& other) const {
    return material == other.material && projectile == other.projectile;
  }
};
struct MaterialHash {
  std::size_t operator()(const MaterialKey& k) const {
    const auto h = std::hash<const void*>{}(k.material);
    return h ^ (std::hash<const void*>{}(k.projectile) + 0x9e3779b9 +
                (h << 6) + (h >> 2));
  }
};
}

ElasticKinematics ComputeElasticKinematics(double m, double M, double T,
                                           double x, double phi) {
  if (!(m > 0 && M > 0 && T >= 0 && x >= 0 && x <= 1) ||
      !std::isfinite(m + M + T + x + phi)) {
    throw std::invalid_argument("Invalid elastic two-body kinematics input");
  }
  ElasticKinematics result;
  const double pLab = std::sqrt(T * (T + 2 * m));
  const double sqrtS = std::sqrt((m + M) * (m + M) + 2 * M * T);
  const double gammaCM = (T + m + M) / sqrtS;
  const double pCM = pLab * M / sqrtS;
  const double transverse = 2 * pCM * std::sqrt(x * (1 - x));
  const double recoilZ = 2 * gammaCM * pCM * x;
  result.momentumCM = pCM;
  result.q2 = 4 * pCM * pCM * x;
  result.recoilKineticEnergy = result.q2 / (2 * M);
  result.finalKineticEnergy = T - result.recoilKineticEnergy;
  if (result.finalKineticEnergy < 0 &&
      result.finalKineticEnergy > -32 * std::numeric_limits<double>::epsilon() * T) {
    result.finalKineticEnergy = 0;
    result.recoilKineticEnergy = T;
  }
  if (result.finalKineticEnergy < 0) {
    throw std::runtime_error("Negative primary kinetic energy after elastic collision");
  }
  result.primaryMomentum = {transverse * std::cos(phi),
                            transverse * std::sin(phi), pLab - recoilZ};
  result.recoilMomentum = {-result.primaryMomentum.x(),
                           -result.primaryMomentum.y(), recoilZ};
  return result;
}

struct G4CCTabulatedIonScatteringModel::Impl {
  struct NucleusState {
    NucleusKey key{};
    double targetMass = 0;
    G4ParticleDefinition* targetDefinition = nullptr;
    const ProjectileTable* table = nullptr;
    double fallbackCrossSectionMass = 0;
    std::unique_ptr<G4IonCoulombCrossSection> fallbackCrossSection;
    std::unique_ptr<G4IonCoulombCrossSection> fallback;
  };
  struct WeightedTarget {
    NucleusState* nucleus = nullptr;
    const G4Element* element = nullptr;
    double numberDensity = 0;
    double sigma = 0;
    EnergyBracket bracket;
  };
  struct MaterialState {
    std::vector<WeightedTarget> targets;
    std::vector<double> cumulative;
    double lastEnergy = -1;
    double total = 0;
  };

  std::shared_ptr<const ElasticTableData> data;
  bool useTables = true;
  double recoilThreshold = 0;
  const G4DataVector* cuts = nullptr;
  G4ParticleChangeForGamma* particleChange = nullptr;
  std::unique_ptr<G4IonCoulombScatteringModel> native;
  std::unordered_map<NucleusKey, std::unique_ptr<NucleusState>, NucleusHash> nuclei;
  std::unordered_map<MaterialKey, std::unique_ptr<MaterialState>, MaterialHash> materials;
  MaterialKey lastMaterialKey{};
  MaterialState* lastMaterialState = nullptr;

  static bool HasCoverage(const MaterialState& material) {
    return std::any_of(material.targets.begin(), material.targets.end(),
                      [](const WeightedTarget& target) {
                        return target.bracket.valid;
                      });
  }

  void InitialiseNative(const G4ParticleDefinition* particle) {
    if (!native) native = std::make_unique<G4IonCoulombScatteringModel>();
    native->SetParticleChange(particleChange);
    native->SetRecoilThreshold(recoilThreshold);
    native->Initialise(particle, *cuts);
  }

  NucleusState& GetNucleus(const G4ParticleDefinition* p, int z, int a, int isomer) {
    const NucleusKey key{p, z, a, isomer};
    auto it = nuclei.find(key);
    if (it != nuclei.end()) return *it->second;
    if (z < 1 || a < z || isomer < 0) {
      Fail("G4CCTabulatedIonScatteringModel::GetNucleus", "Invalid target isotope");
    }
    auto state = std::make_unique<NucleusState>();
    state->key = key;
    state->targetDefinition = G4ParticleTable::GetParticleTable()->GetIonTable()->
        GetIon(z, a, isomer);
    if (!state->targetDefinition) {
      Fail("G4CCTabulatedIonScatteringModel::GetNucleus", "Cannot construct target ion/isomer");
    }
    state->targetMass = state->targetDefinition->GetPDGMass();
    if (useTables && data) {
      state->table = data->FindByZ(p->GetAtomicNumber(), z);
      if(!state->table) WarnMissingPair(p->GetAtomicNumber(), z);
    }
    state->fallbackCrossSectionMass = z == 1
        ? proton_mass_c2 : G4NistManager::Instance()->GetAtomicMassAmu(z) * amu_c2;
    state->fallbackCrossSection = std::make_unique<G4IonCoulombCrossSection>();
    state->fallbackCrossSection->Initialise(p, 1.0);
    state->fallback = std::make_unique<G4IonCoulombCrossSection>();
    state->fallback->Initialise(p, 1.0);
    auto* result = state.get();
    nuclei.emplace(key, std::move(state));
    return *result;
  }

  double CrossSection(NucleusState& state, double energy, EnergyBracket* bracket) {
    if (bracket) *bracket = EnergyBracket{};
    if (!(energy > 0)) return 0;
    if (state.table) {
      const double energyPerNucleon = energy / MeV / state.key.projectile->GetAtomicMass();
      auto b = state.table->Lookup(energyPerNucleon);
      if (b.valid) {
        if (bracket) *bracket = b;
        return b.sigma_mb * millibarn;
      }
    }
    state.fallbackCrossSection->SetupKinematic(
        energy, state.fallbackCrossSectionMass);
    state.fallbackCrossSection->SetupTarget(state.key.z, energy, 0);
    const double sigma = state.fallbackCrossSection->NuclearCrossSection();
    if (!(sigma >= 0) || !std::isfinite(sigma)) {
      Fail("G4CCTabulatedIonScatteringModel::CrossSection", "Invalid stock fallback cross section");
    }
    return sigma;
  }

  MaterialState& GetMaterial(const G4Material* material,
                              const G4ParticleDefinition* p) {
    const MaterialKey key{material, p};
    if (lastMaterialState && key == lastMaterialKey) return *lastMaterialState;
    auto found = materials.find(key);
    if (found != materials.end()) {
      lastMaterialKey = key;
      lastMaterialState = found->second.get();
      return *lastMaterialState;
    }
    auto result = std::make_unique<MaterialState>();
    const auto* atomDensities = material->GetVecNbOfAtomsPerVolume();
    for (std::size_t ei = 0; ei < material->GetNumberOfElements(); ++ei) {
      const auto* element = material->GetElement(ei);
      const auto count = element->GetNumberOfIsotopes();
      const auto* isotopes = element->GetIsotopeVector();
      const auto* abundance = element->GetRelativeAbundanceVector();
      if (!count || !isotopes || !abundance) {
        Fail("G4CCTabulatedIonScatteringModel::GetMaterial",
             "Material element has no explicit/natural isotope composition");
      }
      double sum = 0;
      for (std::size_t ii = 0; ii < count; ++ii) sum += abundance[ii];
      if (!(sum > 0) || !std::isfinite(sum)) {
        Fail("G4CCTabulatedIonScatteringModel::GetMaterial", "Invalid isotope abundances");
      }
      for (std::size_t ii = 0; ii < count; ++ii) {
        if (abundance[ii] == 0) continue;
        if (!(abundance[ii] >= 0) || !std::isfinite(abundance[ii])) {
          Fail("G4CCTabulatedIonScatteringModel::GetMaterial", "Invalid isotope abundance");
        }
        const auto* isotope = (*isotopes)[ii];
        auto& state = GetNucleus(p, isotope->GetZ(), isotope->GetN(), isotope->Getm());
        result->targets.push_back({&state, element,
                                  atomDensities[ei] * abundance[ii] / sum, 0, {}});
      }
    }
    result->cumulative.resize(result->targets.size());
    lastMaterialKey = key;
    lastMaterialState = result.get();
    materials.emplace(key, std::move(result));
    return *lastMaterialState;
  }

  double UpdateWeights(MaterialState& material, double energy) {
    if (material.lastEnergy == energy) return material.total;
    material.total = 0;
    for (std::size_t i = 0; i < material.targets.size(); ++i) {
      auto& t = material.targets[i];
      t.sigma = CrossSection(*t.nucleus, energy, &t.bracket);
      material.total += t.numberDensity * t.sigma;
      material.cumulative[i] = material.total;
    }
    material.lastEnergy = energy;
    return material.total;
  }
};

G4CCTabulatedIonScatteringModel::G4CCTabulatedIonScatteringModel(
    std::shared_ptr<const ElasticTableData> data, bool useTables)
    : G4VEmModel(useTables ? "CCTabulatedIonScatteringModel" : "CCNativeIonCoulombFormula"),
      impl_(std::make_unique<Impl>()) {
  impl_->data = std::move(data);
  impl_->useTables = useTables;
}

G4CCTabulatedIonScatteringModel::~G4CCTabulatedIonScatteringModel() = default;

void G4CCTabulatedIonScatteringModel::Initialise(const G4ParticleDefinition* particle,
                                         const G4DataVector& cuts) {
  impl_->cuts = &cuts;
  impl_->particleChange = GetParticleChangeForGamma();
  impl_->InitialiseNative(particle);
  impl_->materials.clear();
  impl_->nuclei.clear();
  impl_->lastMaterialState = nullptr;
}

void G4CCTabulatedIonScatteringModel::InitialiseLocal(const G4ParticleDefinition* particle,
                                              G4VEmModel* master) {
  auto* model = dynamic_cast<G4CCTabulatedIonScatteringModel*>(master);
  if (!model) Fail("G4CCTabulatedIonScatteringModel::InitialiseLocal", "Incorrect master model");
  impl_->data = model->impl_->data;
  impl_->useTables = model->impl_->useTables;
  impl_->recoilThreshold = model->impl_->recoilThreshold;
  if (!impl_->cuts) impl_->cuts = model->impl_->cuts;
  impl_->particleChange = GetParticleChangeForGamma();
  impl_->InitialiseNative(particle);
}

void G4CCTabulatedIonScatteringModel::StartTracking(G4Track* track) {
  if (impl_->native) impl_->native->StartTracking(track);
}

void G4CCTabulatedIonScatteringModel::SetRecoilThreshold(G4double value) {
  if (!(value >= 0) || !std::isfinite(value)) {
    Fail("G4CCTabulatedIonScatteringModel::SetRecoilThreshold", "Invalid recoil threshold");
  }
  impl_->recoilThreshold = value;
  if (impl_->native) impl_->native->SetRecoilThreshold(value);
}

G4double G4CCTabulatedIonScatteringModel::CrossSectionPerVolume(
    const G4Material* material, const G4ParticleDefinition* p, G4double energy,
    G4double cut, G4double maximum) {
  try {
    auto& state = impl_->GetMaterial(material, p);
    const double total = impl_->UpdateWeights(state, energy);
    if (!Impl::HasCoverage(state) && CurrentCouple() && impl_->native) {
      impl_->native->SetCurrentCouple(CurrentCouple());
      return impl_->native->CrossSectionPerVolume(material, p, energy,
                                                 cut, maximum);
    }
    return total;
  } catch (const std::exception& e) {
    Fail("G4CCTabulatedIonScatteringModel::CrossSectionPerVolume", e.what());
  }
}

G4double G4CCTabulatedIonScatteringModel::CrossSectionForIsotope(
    const G4ParticleDefinition* p, G4double energy, int z, int a, int isomer) {
  try {
    auto& state = impl_->GetNucleus(p, z, a, isomer);
    return impl_->CrossSection(state, energy, nullptr);
  } catch (const std::exception& e) {
    Fail("G4CCTabulatedIonScatteringModel::CrossSectionForIsotope", e.what());
  }
}

G4double G4CCTabulatedIonScatteringModel::ComputeCrossSectionPerAtom(
    const G4ParticleDefinition* p, G4double energy, G4double Z, G4double A,
    G4double, G4double) {
  const int z = static_cast<int>(std::lround(Z));
  const int a = static_cast<int>(std::lround(A));
  if (a >= z && std::abs(A - a) <= 1e-8) {
    return CrossSectionForIsotope(p, energy, z, a);
  }
  auto* element = G4NistManager::Instance()->FindOrBuildElement(z, true);
  if (!element) Fail("G4CCTabulatedIonScatteringModel::ComputeCrossSectionPerAtom", "Unknown element");
  const auto* isotopes = element->GetIsotopeVector();
  const auto* abundances = element->GetRelativeAbundanceVector();
  double sum = 0, sigma = 0;
  for (std::size_t i = 0; i < element->GetNumberOfIsotopes(); ++i) {
    const auto* isotope = (*isotopes)[i];
    sum += abundances[i];
    sigma += abundances[i] * CrossSectionForIsotope(
        p, energy, isotope->GetZ(), isotope->GetN(), isotope->Getm());
  }
  if (!(sum > 0)) Fail("G4CCTabulatedIonScatteringModel::ComputeCrossSectionPerAtom", "Empty isotope mixture");
  return sigma / sum;
}

void G4CCTabulatedIonScatteringModel::SampleSecondaries(
    std::vector<G4DynamicParticle*>* secondaries,
    const G4MaterialCutsCouple* couple, const G4DynamicParticle* dp,
    G4double minimum, G4double maximum) {
  try {
    const double energy = dp->GetKineticEnergy();
    const auto* p = dp->GetDefinition();
    auto& material = impl_->GetMaterial(couple->GetMaterial(), p);
    const double total = impl_->UpdateWeights(material, energy);
    if (!(total > 0)) return;
    if (!Impl::HasCoverage(material)) {
      if (!impl_->native || !impl_->particleChange) {
        Fail("G4CCTabulatedIonScatteringModel::SampleSecondaries",
             "Native delegate is not initialised");
      }
      impl_->native->SetCurrentCouple(couple);
      impl_->native->SampleSecondaries(secondaries, couple, dp, minimum, maximum);
      const auto* element = impl_->native->GetCurrentElement();
      SetCurrentElement(element);
      return;
    }
    const double targetUniform = G4UniformRand() * total;
    auto selected = std::upper_bound(material.cumulative.begin(),
                                     material.cumulative.end(), targetUniform);
    const auto index = std::min<std::size_t>(
        selected - material.cumulative.begin(), material.targets.size() - 1);
    auto& target = material.targets[index];
    auto& nucleus = *target.nucleus;
    SetCurrentElement(target.element);

    double x;
    int node = -1;
    if (target.bracket.valid) {
      node = static_cast<int>(nucleus.table->SelectNode(target.bracket, G4UniformRand()));
      x = nucleus.table->SampleNodeX(static_cast<std::size_t>(node), G4UniformRand());
    } else {
      nucleus.fallback->SetupKinematic(energy, nucleus.targetMass);
      nucleus.fallback->SetupTarget(nucleus.key.z, energy, 0);
      x = 0.5 * std::clamp(nucleus.fallback->SampleCosineTheta(), 0.0, 2.0);
    }
    const double phi = twopi * G4UniformRand();
    const double projectileMass = dp->GetMass();
    const auto k = ComputeElasticKinematics(projectileMass, nucleus.targetMass,
                                            energy, x, phi);
    const auto incomingDirection = dp->GetMomentumDirection();
    auto newDirection = k.primaryMomentum.mag2() > 0
        ? k.primaryMomentum.unit() : G4ThreeVector(0, 0, 1);
    newDirection.rotateUz(incomingDirection);
    auto* change = impl_->particleChange;
    if (!change) Fail("G4CCTabulatedIonScatteringModel::SampleSecondaries", "Model is not initialised");
    change->ProposeMomentumDirection(newDirection);
    change->SetProposedKineticEnergy(k.finalKineticEnergy);

    double cut = impl_->recoilThreshold;
    if (impl_->cuts) {
      const auto ci = static_cast<std::size_t>(couple->GetIndex());
      if (ci >= impl_->cuts->size()) {
        Fail("G4CCTabulatedIonScatteringModel::SampleSecondaries", "Missing recoil production cut");
      }
      cut = std::max(cut, (*impl_->cuts)[ci]);
    }
    const bool recoilCreated = k.recoilKineticEnergy > cut;
    if (recoilCreated) {
      auto recoilDirection = k.recoilMomentum.unit();
      recoilDirection.rotateUz(incomingDirection);
      secondaries->push_back(new G4DynamicParticle(nucleus.targetDefinition,
                                                   recoilDirection,
                                                   k.recoilKineticEnergy));
    }
    const double localDeposit = recoilCreated ? 0 : k.recoilKineticEnergy;
    change->ProposeLocalEnergyDeposit(localDeposit);
    change->ProposeNonIonizingEnergyDeposit(localDeposit);

  } catch (const std::exception& e) {
    Fail("G4CCTabulatedIonScatteringModel::SampleSecondaries", e.what());
  }
}

void G4CCTabulatedIonScatteringModel::ModelDescription(std::ostream& out) const {
  out << "Full-range elastic scattering with tables selected by projectile Z "
         "and target Z for every isotope and internal state, at actual T/A. "
         "Tabulated screened "
         "sigma and exact within-cell R(x)/(As+x)^2 sampling; between-energy "
         "interpolation is a positive linear mixture of unnormalised measures. "
         "Wholly uncovered materials/energies delegate to the untouched "
         "G4IonCoulombScatteringModel for sigma and final state. Mixed "
         "covered/uncovered materials select the target only once; uncovered "
         "targets then use G4IonCoulombCrossSection, NIST mean mass for sigma "
         "and actual isotope mass for the angular sampler. Primary charge "
         "and internal state are preserved.\n";
}
}
