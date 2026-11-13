#include "G4CCAtimaEquilibriumDedxTable.hh"
#include "G4CCCurrentChargeConfig.hh"

#include "G4CCIonChargeState.hh"
#include "G4Element.hh"
#include "G4Exception.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
  constexpr G4double kMinEnergy = 2.0;
  constexpr G4double kMaxEnergy = 40.0;
  constexpr std::size_t kExpectedCurves = 68*93;

  [[noreturn]] void InvalidTable(const G4String& path, const G4String& reason)
  {
    const G4String message = "Invalid equilibrium ATIMA table '" + path + "': " + reason;
    G4Exception("G4CCAtimaEquilibriumDedxTable", "CCAtimaEquilibrium001",
                FatalException, message.c_str());
    throw std::runtime_error(message);
  }

  G4String Trim(const G4String& text)
  {
    std::size_t begin = 0, end = text.size();
    while(begin < end && std::isspace(static_cast<unsigned char>(text[begin]))) { ++begin; }
    while(end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) { --end; }
    return text.substr(begin, end - begin);
  }

  template<std::size_t N>
  G4bool Split(const G4String& line, std::array<G4String, N>& fields)
  {
    std::size_t begin = 0;
    for(std::size_t i = 0; i < fields.size(); ++i) {
      const auto end = line.find(',', begin);
      if((i + 1 < fields.size()) != (end != G4String::npos)) { return false; }
      fields[i] = Trim(line.substr(begin, end == G4String::npos ? end : end - begin));
      if(fields[i].empty()) { return false; }
      begin = end == G4String::npos ? line.size() : end + 1;
    }
    return true;
  }

  G4bool Integer(const G4String& text, G4int& value)
  {
    try {
      std::size_t used = 0;
      const long parsed = std::stol(text, &used);
      if(used != text.size() || parsed < 1 || parsed > 96) { return false; }
      value = static_cast<G4int>(parsed);
      return true;
    } catch(const std::exception&) { return false; }
  }

  G4bool PositiveDouble(const G4String& text, G4double& value)
  {
    try {
      std::size_t used = 0;
      value = std::stod(text, &used);
      return used == text.size() && std::isfinite(value) && value > 0.0;
    } catch(const std::exception&) { return false; }
  }

  G4bool SupportedTarget(G4int zt)
  {
    return zt >= 1 && zt <= 96 && !(zt >= 85 && zt <= 87);
  }
}

const G4CCAtimaEquilibriumDedxTable& G4CCAtimaEquilibriumDedxTable::Instance()
{
  static const G4CCAtimaEquilibriumDedxTable table;
  return table;
}

G4CCAtimaEquilibriumDedxTable::G4CCAtimaEquilibriumDedxTable()
  : fPath(SelectPath()), fData(Load(fPath))
{}

G4String G4CCAtimaEquilibriumDedxTable::SelectPath()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode())
    return G4CCCurrentChargeConfig::Instance().Options().tables.atimaLow;
  const char* path = std::getenv("G4CC_ATIMA_EQUILIBRIUM_DEDX_TABLE");
  if(path != nullptr && path[0] != '\0') { return path; }
  InvalidTable("", "Set G4CC_ATIMA_EQUILIBRIUM_DEDX_TABLE to the equilibrium ATIMA CSV "
    "for Legacy mode (for example, runtime-tables/atima_low.csv).");
}

G4CCAtimaEquilibriumDedxTable::Data
G4CCAtimaEquilibriumDedxTable::Load(const G4String& path)
{
  std::ifstream input(path.c_str());
  if(!input) { InvalidTable(path, "cannot open selected CSV"); }
  G4String line;
  if(!std::getline(input, line)) { InvalidTable(path, "missing CSV header"); }
  if(line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xef &&
     static_cast<unsigned char>(line[1]) == 0xbb &&
     static_cast<unsigned char>(line[2]) == 0xbf) { line.erase(0, 3); }
  std::array<G4String, 5> legacyFields;
  std::array<G4String, 6> fields;
  const std::array<G4String, 5> expected = {{"Zp", "Zt", "E_MeV_u",
    "dEdx_MeV_cm2_g", "dEdx_left_MeV_cm2_g"}};
  const std::array<G4String, 6> expectedReference = {{"Zp", "Zt", "At_ref_g_mol", "E_MeV_u",
    "dEdx_MeV_cm2_g", "dEdx_left_MeV_cm2_g"}};
  const bool referenceFormat = Split(line, fields) && fields == expectedReference;
  const bool legacyFormat = Split(line, legacyFields) && legacyFields == expected;
  const bool explicitMode = G4CCCurrentChargeConfig::Instance().IsDiracMode();
  if(!referenceFormat && (!legacyFormat || explicitMode)) {
    InvalidTable(path, "explicit Dirac mode requires Zp,Zt,At_ref_g_mol,E_MeV_u,"
      "dEdx_MeV_cm2_g,dEdx_left_MeV_cm2_g. Upgrade legacy CSV using convert_low_atima_reference.py.");
  }

  Data data;
  std::size_t lineNumber = 1;
  G4int previousKey = -1;
  G4double previousEnergy = 0.0;
  while(std::getline(input, line)) {
    ++lineNumber;
    if(Trim(line).empty()) { continue; }
    G4int zp = 0, zt = 0;
    G4double referenceMass = 0.0;
    Point point{};
    const bool valid = referenceFormat
      ? (Split(line, fields) && Integer(fields[0], zp) && Integer(fields[1], zt) &&
         PositiveDouble(fields[2], referenceMass) && PositiveDouble(fields[3], point.energy) &&
         PositiveDouble(fields[4], point.main) && PositiveDouble(fields[5], point.left))
      : (Split(line, legacyFields) && Integer(legacyFields[0], zp) && Integer(legacyFields[1], zt) &&
         PositiveDouble(legacyFields[2], point.energy) && PositiveDouble(legacyFields[3], point.main) &&
         PositiveDouble(legacyFields[4], point.left));
    if(!valid) {
      InvalidTable(path, "malformed/nonpositive/nonfinite row at line " + std::to_string(lineNumber));
    }
    if(zp < 29 || !SupportedTarget(zt) ||
       point.energy < kMinEnergy || point.energy > kMaxEnergy) {
      InvalidTable(path, "row outside the declared Z/E domain at line " + std::to_string(lineNumber));
    }
    const G4int key = Key(zp, zt);
    if(key < previousKey || (key == previousKey && point.energy <= previousEnergy)) {
      InvalidTable(path, "duplicate or unsorted (Zp,Zt,E) at line " + std::to_string(lineNumber));
    }
    auto& curve = data[key];
    if(curve.points.empty()) { curve.referenceAtomicMass = referenceMass; }
    else if(curve.referenceAtomicMass != referenceMass) {
      InvalidTable(path, "inconsistent reference target mass at line " + std::to_string(lineNumber));
    }
    curve.points.push_back(point);
    previousKey = key;
    previousEnergy = point.energy;
  }
  if(input.bad()) { InvalidTable(path, "I/O error while reading CSV"); }
  if(data.empty()) { InvalidTable(path, "no projectile-target curves"); }
  if(!explicitMode && data.size() != kExpectedCurves) {
    InvalidTable(path, "legacy mode expected all 6324 projectile-target curves");
  }
  for(const auto& entry : data) {
    const auto& points = entry.second.points;
    if(points.size() < 2 || points.front().energy != kMinEnergy ||
       points.back().energy != kMaxEnergy) {
      InvalidTable(path, "incomplete 2--40 curve for Zp=" + std::to_string(entry.first/1000) +
                         ", Zt=" + std::to_string(entry.first%1000));
    }
  }
  return data;
}

std::size_t G4CCAtimaEquilibriumDedxTable::NumberOfRows() const
{
  std::size_t count = 0;
  for(const auto& entry : fData) { count += entry.second.points.size(); }
  return count;
}

G4bool G4CCAtimaEquilibriumDedxTable::Lookup(G4int zp, G4int zt, G4double energy,
                                            G4bool leftLimit, G4double& value) const
{
  value = 0.0;
  if(!std::isfinite(energy) || energy < kMinEnergy || energy > kMaxEnergy ||
     zp < 29 || zp > 96) { return false; }
  if(!SupportedTarget(zt)) {
    InvalidTable(fPath, "unsupported selected target Zt=" + std::to_string(zt));
  }
  const auto found = fData.find(Key(zp, zt));
  if(found == fData.end()) {
    InvalidTable(fPath, "missing selected curve for Zp=" + std::to_string(zp) +
                       ", Zt=" + std::to_string(zt));
  }
  const auto& curve = found->second.points;
  const auto upper = std::lower_bound(curve.begin(), curve.end(), energy,
    [](const Point& point, G4double e) { return point.energy < e; });
  if(upper != curve.end() && upper->energy == energy) {
    value = leftLimit ? upper->left : upper->main;
    return true;
  }
  if(upper == curve.begin() || upper == curve.end()) {
    InvalidTable(fPath, "selected energy is not bracketed by the validated curve");
  }
  const Point& lo = *(upper - 1);
  const G4double x = std::log1p((energy - lo.energy)/lo.energy) /
                     std::log1p((upper->energy - lo.energy)/lo.energy);
  value = std::exp((1.0 - x)*std::log(lo.main) + x*std::log(upper->left));
  if(!std::isfinite(value) || !(value > 0.0)) {
    InvalidTable(fPath, "nonfinite/nonpositive interpolated stopping");
  }
  return true;
}

G4bool G4CCAtimaEquilibriumDedxTable::GetMassDEDX(G4int zp, G4int zt,
                                               G4double energy, G4double& value) const
{
  return Lookup(zp, zt, energy, false, value);
}

G4bool G4CCAtimaEquilibriumDedxTable::HasPair(G4int zp, G4int zt) const
{
  return fData.find(Key(zp,zt)) != fData.end();
}

void G4CCAtimaEquilibriumDedxTable::RequirePair(G4int zp, G4int zt) const
{
  if(!HasPair(zp,zt)) {
    InvalidTable(fPath, "missing required low-energy curve for Zp=" + std::to_string(zp) +
                       ", Zt=" + std::to_string(zt));
  }
}

G4double G4CCAtimaEquilibriumDedxTable::GetReferenceAtomicMass(G4int zp, G4int zt) const
{
  RequirePair(zp,zt);
  const double mass = fData.at(Key(zp,zt)).referenceAtomicMass;
  if(!(mass > 0.0)) { InvalidTable(fPath, "legacy curve has no reference target atomic mass"); }
  return mass;
}

G4bool G4CCAtimaEquilibriumDedxTable::GetMassDEDXLeftLimit(G4int zp, G4int zt,
                                                        G4double energy, G4double& value) const
{
  return Lookup(zp, zt, energy, true, value);
}

G4bool G4CCAtimaEquilibriumDedxTable::GetMaterialDEDXPerVolume(
  const G4ParticleDefinition* particle, const G4Material* material,
  G4double kineticEnergy, G4double& dedxPerVolume) const
{
  dedxPerVolume = 0.0;
  if(particle == nullptr || material == nullptr) { return false; }
  const G4int zp = G4CCIonChargeState::GetZp(particle);
  const G4double energy = G4CCIonChargeState::KineticEnergyPerNucleon(particle, kineticEnergy)/CLHEP::MeV;
  if(!std::isfinite(energy) || energy < kMinEnergy || energy > kMaxEnergy ||
     zp < 29 || zp > 96) { return false; }
  const auto* elements = material->GetElementVector();
  const G4double* fractions = material->GetFractionVector();
  const G4double* atomDensities = material->GetVecNbOfAtomsPerVolume();
  const auto count = material->GetNumberOfElements();
  if(elements == nullptr || fractions == nullptr || atomDensities == nullptr || count == 0) {
    InvalidTable(fPath, "selected material has no element densities/fractions");
  }
  G4double massStopping = 0.0;
  for(std::size_t i = 0; i < count; ++i) {
    if(!std::isfinite(fractions[i]) || fractions[i] < 0.0) {
      InvalidTable(fPath, "invalid selected material mass fraction");
    }
    if(fractions[i] == 0.0) { continue; }
    G4double elemental = 0.0;
    const G4int zt = static_cast<G4int>(std::lround((*elements)[i]->GetZ()));
    GetMassDEDX(zp, zt, energy, elemental);
    const double referenceMass = fData.at(Key(zp,zt)).referenceAtomicMass;
    if(referenceMass > 0.0) {
      if(!std::isfinite(atomDensities[i]) || atomDensities[i] < 0.0) {
        InvalidTable(fPath, "invalid selected material atomic density");
      }
      dedxPerVolume += atomDensities[i]*elemental*referenceMass /
        (CLHEP::Avogadro*CLHEP::mole)*(CLHEP::MeV*CLHEP::cm2);
    } else { massStopping += fractions[i]*elemental; }
  }
  dedxPerVolume += massStopping*(CLHEP::MeV*CLHEP::cm2/CLHEP::g)*material->GetDensity();
  if(!std::isfinite(dedxPerVolume) || !(dedxPerVolume > 0.0)) {
    InvalidTable(fPath, "nonpositive/nonfinite selected material stopping");
  }
  return true;
}
