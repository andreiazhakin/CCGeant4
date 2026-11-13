#include "G4CCAtimaQDedxTable.hh"

#include "G4CCAtimaEquilibriumDedxTable.hh"
#include "G4AutoLock.hh"
#include "G4Exception.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4CCIonChargeState.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

namespace
{
  G4Mutex g4ccAtimaQDedxMutex = G4MUTEX_INITIALIZER;

  G4String TrimCopy(const G4String& in)
  {
    std::size_t b = 0;
    while(b < in.size() && std::isspace(static_cast<unsigned char>(in[b]))) { ++b; }
    std::size_t e = in.size();
    while(e > b && std::isspace(static_cast<unsigned char>(in[e - 1]))) { --e; }
    return in.substr(b, e - b);
  }

  std::vector<G4String> SplitCSVLine(const G4String& line)
  {
    std::vector<G4String> out;
    G4String current;
    bool quoted = false;
    for(char c : line) {
      if(c == '"') {
        quoted = !quoted;
      } else if(c == ',' && !quoted) {
        out.push_back(TrimCopy(current));
        current.clear();
      } else {
        current.push_back(c);
      }
    }
    out.push_back(TrimCopy(current));
    return out;
  }

  G4int FindColumn(const std::vector<G4String>& header, const std::vector<G4String>& names)
  {
    for(const auto& name : names) {
      for(std::size_t i = 0; i < header.size(); ++i) {
        if(header[i] == name) { return static_cast<G4int>(i); }
      }
    }
    return -1;
  }

  G4bool ToDouble(const std::vector<G4String>& fields, G4int idx, G4double& value)
  {
    if(idx < 0 || static_cast<std::size_t>(idx) >= fields.size()) { return false; }
    try {
      value = std::stod(fields[static_cast<std::size_t>(idx)]);
    } catch(...) {
      return false;
    }
    return std::isfinite(value);
  }

  G4bool ToInt(const std::vector<G4String>& fields, G4int idx, G4int& value)
  {
    G4double v = 0.0;
    if(!ToDouble(fields, idx, v)) { return false; }
    value = static_cast<G4int>(std::lround(v));
    return true;
  }
}

G4CCAtimaQDedxTable& G4CCAtimaQDedxTable::Instance()
{
  static G4CCAtimaQDedxTable table;
  return table;
}

G4CCAtimaQDedxTable::G4CCAtimaQDedxTable()
{
  const char* path = std::getenv("G4CC_ATIMA_Q_DEDX_TABLE");
  if(path == nullptr) { path = std::getenv("G4CC_ATIMA_DEDX_TABLE"); }
  if(path != nullptr) { fPath = path; }
}

void G4CCAtimaQDedxTable::CheckConfigurationChangeAllowed()
{
  const auto state = G4StateManager::GetStateManager()->GetCurrentState();
  if(G4Threading::IsWorkerThread() ||
     (state != G4State_PreInit && state != G4State_Idle)) {
    const char* message = "ATIMA Q table configuration may change only on the "
      "master in PreInit/Idle, with no table readers active.";
    G4Exception("G4CCAtimaQDedxTable::CheckConfigurationChangeAllowed",
                "CCAtimaQConfig001", FatalException, message);
    throw std::runtime_error(message);
  }
}

const G4String& G4CCAtimaQDedxTable::Path() const
{
  return fPath;
}

void G4CCAtimaQDedxTable::SetPath(const G4String& path)
{
  CheckConfigurationChangeAllowed();
  G4AutoLock lock(&g4ccAtimaQDedxMutex);
  fPath = path;
  Clear();
}

G4bool G4CCAtimaQDedxTable::IsLoaded() const
{
  return fLoaded.load(std::memory_order_acquire);
}

G4bool G4CCAtimaQDedxTable::IsAvailable() const
{
  return EnsureLoaded();
}

void G4CCAtimaQDedxTable::Clear()
{
  fData.clear();
  fLoaded.store(false, std::memory_order_release);
  fLoadAttempted.store(false, std::memory_order_release);
}

long G4CCAtimaQDedxTable::Key(G4int zProjectile, G4int zTarget)
{
  return 1000L*static_cast<long>(zProjectile) + static_cast<long>(zTarget);
}

G4bool G4CCAtimaQDedxTable::EnsureLoaded() const
{
  if(fLoadAttempted.load(std::memory_order_acquire)) {
    return fLoaded.load(std::memory_order_acquire);
  }

  G4AutoLock lock(&g4ccAtimaQDedxMutex);
  if(fPath.empty()) { return false; }
  if(!fLoadAttempted.load(std::memory_order_acquire)) {
    const G4bool loaded = Load();
    fLoaded.store(loaded, std::memory_order_release);
    fLoadAttempted.store(true, std::memory_order_release);
  }
  return fLoaded.load(std::memory_order_acquire);
}

G4bool G4CCAtimaQDedxTable::Load() const
{
  fData.clear();

  std::ifstream input(fPath.c_str());
  if(!input) {
    G4cout << "G4CCAtimaQDedxTable: cannot open table: " << fPath << G4endl;
    return false;
  }

  G4String line;
  if(!std::getline(input, line)) { return false; }
  auto header = SplitCSVLine(line);

  const G4int cZp = FindColumn(header, {"Zp", "zProjectile"});
  const G4int cZt = FindColumn(header, {"Zt", "zTarget"});
  const G4int cE = FindColumn(header, {"E_MeV_u", "E_MeV_per_u", "energy_MeV_u"});
  const G4int cEqDedx = FindColumn(header, {"ATIMA_dEdx_eq_MeV_cm2_g",
                                             "ATIMA_full_dEdx_MeV_cm2_g",
                                             "ATIMA_full_MeV_cm2_g"});
  const G4int cQ2Eq = FindColumn(header, {"q2_eq", "q_rms2_eq"});
  const G4int cQrms = FindColumn(header, {"q_rms_eq", "qrms_eq"});

  const G4int cQ = FindColumn(header, {"q", "charge"});
  const G4int cDedxQ = FindColumn(header, {"dedx_q_MeV_cm2_g",
                                           "dEdx_q_MeV_cm2_g"});

  const G4bool compactFormat = (cZp >= 0 && cZt >= 0 && cE >= 0 && cEqDedx >= 0 &&
                                (cQ2Eq >= 0 || cQrms >= 0));
  const G4bool qResolvedFormat = (cZp >= 0 && cZt >= 0 && cE >= 0 && cQ >= 0 && cDedxQ >= 0);

  if(!compactFormat && !qResolvedFormat) {
    G4cout << "G4CCAtimaQDedxTable: unsupported CSV columns in " << fPath << G4endl;
    return false;
  }

  std::size_t nRows = 0;
  while(std::getline(input, line)) {
    if(line.empty() || line[0] == '#') { continue; }
    auto fields = SplitCSVLine(line);

    G4int zp = 0;
    G4int zt = 0;
    G4double eMeVPerU = 0.0;
    if(!ToInt(fields, cZp, zp) || !ToInt(fields, cZt, zt) || !ToDouble(fields, cE, eMeVPerU)) {
      continue;
    }

    G4double unitDedx = 0.0;
    if(compactFormat) {
      G4double eqDedx = 0.0;
      G4double q2eq = 0.0;
      if(!ToDouble(fields, cEqDedx, eqDedx)) { continue; }
      if(cQ2Eq >= 0) {
        if(!ToDouble(fields, cQ2Eq, q2eq)) { continue; }
      } else {
        G4double qrms = 0.0;
        if(!ToDouble(fields, cQrms, qrms)) { continue; }
        q2eq = qrms*qrms;
      }
      if(!(eqDedx > 0.0) || !(q2eq > 0.0)) { continue; }
      unitDedx = eqDedx/q2eq;
    } else {
      G4int q = 0;
      if(!ToInt(fields, cQ, q) || q != 1) { continue; }
      if(!ToDouble(fields, cDedxQ, unitDedx) || !(unitDedx > 0.0)) { continue; }
    }

    fData[Key(zp, zt)].push_back({eMeVPerU, unitDedx});
    ++nRows;
  }

  for(auto& item : fData) {
    auto& curve = item.second;
    std::sort(curve.begin(), curve.end(), [](const Point& a, const Point& b) {
      return a.eMeVPerU < b.eMeVPerU;
    });
    curve.erase(std::unique(curve.begin(), curve.end(), [](const Point& a, const Point& b) {
      return a.eMeVPerU == b.eMeVPerU;
    }), curve.end());
  }

  if(fData.empty()) {
    G4cout << "G4CCAtimaQDedxTable: no usable rows in " << fPath << G4endl;
    return false;
  }

  G4cout << "G4CCAtimaQDedxTable: loaded " << nRows
         << " rows, " << fData.size() << " Zp-Zt curves from " << fPath << G4endl;
  return true;
}

G4bool G4CCAtimaQDedxTable::GetUnitMassDEDX(G4int zProjectile,
                                            G4int zTarget,
                                            G4double kineticEnergyMeVPerU,
                                            G4double& valueMeVcm2PerG) const
{
  valueMeVcm2PerG = 0.0;
  if(!EnsureLoaded() || kineticEnergyMeVPerU <= 0.0) { return false; }

  const auto it = fData.find(Key(zProjectile, zTarget));
  if(it == fData.end() || it->second.empty()) { return false; }
  const Curve& curve = it->second;

  if(kineticEnergyMeVPerU < curve.front().eMeVPerU ||
     kineticEnergyMeVPerU > curve.back().eMeVPerU) {
    return false;
  }

  if(curve.size() == 1 || kineticEnergyMeVPerU == curve.front().eMeVPerU) {
    valueMeVcm2PerG = curve.front().unitMassDedx;
    return valueMeVcm2PerG > 0.0;
  }
  if(kineticEnergyMeVPerU == curve.back().eMeVPerU) {
    valueMeVcm2PerG = curve.back().unitMassDedx;
    return valueMeVcm2PerG > 0.0;
  }

  auto upper = std::upper_bound(curve.begin(), curve.end(), kineticEnergyMeVPerU,
                                [](G4double e, const Point& p) { return e < p.eMeVPerU; });
  if(upper == curve.begin() || upper == curve.end()) { return false; }
  const Point& p1 = *upper;
  const Point& p0 = *(upper - 1);
  if(!(p0.eMeVPerU > 0.0) || !(p1.eMeVPerU > p0.eMeVPerU) ||
     !(p0.unitMassDedx > 0.0) || !(p1.unitMassDedx > 0.0)) {
    return false;
  }

  const G4double x = (std::log(kineticEnergyMeVPerU) - std::log(p0.eMeVPerU)) /
                     (std::log(p1.eMeVPerU) - std::log(p0.eMeVPerU));
  valueMeVcm2PerG = std::exp(std::log(p0.unitMassDedx) + x*(std::log(p1.unitMassDedx) - std::log(p0.unitMassDedx)));
  return valueMeVcm2PerG > 0.0;
}

G4bool G4CCAtimaQDedxTable::GetMaterialUnitMassDEDX(const G4ParticleDefinition* particle,
                                                    const G4Material* material,
                                                    G4double kineticEnergy,
                                                    G4double& valueMeVcm2PerG) const
{
  valueMeVcm2PerG = 0.0;
  if(particle == nullptr || material == nullptr) { return false; }

  const G4int zp = G4CCIonChargeState::GetZp(particle);
  if(zp <= 0) { return false; }

  const G4double eMeVPerU = G4CCIonChargeState::KineticEnergyPerNucleon(particle, kineticEnergy)/CLHEP::MeV;
  if(eMeVPerU <= 0.0) { return false; }

  const auto* elements = material->GetElementVector();
  const G4double* fractions = material->GetFractionVector();
  const G4int n = material->GetNumberOfElements();
  if(elements == nullptr || fractions == nullptr || n <= 0) { return false; }

  G4double sum = 0.0;
  for(G4int i = 0; i < n; ++i) {
    const G4int zt = static_cast<G4int>(std::lround((*elements)[i]->GetZ()));
    G4double elemental = 0.0;
    if(!GetUnitMassDEDX(zp, zt, eMeVPerU, elemental)) { return false; }
    const G4double referenceMass =
      G4CCAtimaEquilibriumDedxTable::Instance().GetReferenceAtomicMass(zp, zt);
    const G4double actualMass = (*elements)[i]->GetA()/(CLHEP::g/CLHEP::mole);
    sum += fractions[i]*elemental*referenceMass/actualMass;
  }

  valueMeVcm2PerG = sum;
  return valueMeVcm2PerG > 0.0;
}

G4bool G4CCAtimaQDedxTable::GetMaterialDEDXPerVolume(const G4ParticleDefinition* particle,
                                                     const G4Material* material,
                                                     G4double kineticEnergy,
                                                     G4double chargeSquare,
                                                     G4double& dedxPerVolume) const
{
  dedxPerVolume = 0.0;
  G4double unitMassDedx = 0.0;
  if(!GetMaterialUnitMassDEDX(particle, material, kineticEnergy, unitMassDedx)) { return false; }
  if(chargeSquare < 0.0) { chargeSquare = 0.0; }

  dedxPerVolume = unitMassDedx*chargeSquare*(CLHEP::MeV*CLHEP::cm2/CLHEP::g)*material->GetDensity();
  return dedxPerVolume > 0.0;
}
