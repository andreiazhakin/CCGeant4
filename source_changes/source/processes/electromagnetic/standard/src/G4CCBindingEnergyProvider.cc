#include "G4CCBindingEnergyProvider.hh"

#include "G4CCCurrentChargeConfig.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>

namespace {

G4double BindingOverrideEnergy()
{
  static const G4double value = []() {
    const char* envMeV = std::getenv("G4CC_BINDING_OVERRIDE_MEV");
    if(envMeV != nullptr) {
      const G4double v = std::atof(envMeV);
      if(v > 0.0) { return v*MeV; }
    }
    const char* envGeV = std::getenv("G4CC_BINDING_OVERRIDE_GEV");
    if(envGeV != nullptr) {
      const G4double v = std::atof(envGeV);
      if(v > 0.0) { return v*GeV; }
    }
    if(envMeV == nullptr && envGeV == nullptr) { return -1.0; }
    return 0.0;
  }();
  return value < 0.0
    ? (G4CCCurrentChargeConfig::Instance().IsDiracMode() ? 0.0 : 1000.0*GeV)
    : value;
}

struct ShellSpec { G4int n; G4int capacity; G4double relFactor; };

constexpr std::array<ShellSpec, 12> kShells{{
  {1, 2, 1.13},
  {2, 2, 1.02},
  {2, 6, 0.96},
  {3, 2, 1.00},
  {3, 6, 0.96},
  {3,10, 0.92},
  {4, 2, 1.00},
  {4, 6, 0.96},
  {4,10, 0.92},
  {4,14, 0.88},
  {5, 2, 1.00},
  {5, 6, 0.96}
}};
}

G4double G4CCHydrogenicBindingEnergyProvider::BindingEnergy(G4int Zp, G4int electronIndex) const
{
  if(Zp <= 0 || electronIndex < 0) { return 0.0; }

  const G4double overrideEnergy = BindingOverrideEnergy();
  if(overrideEnergy > 0.0) { return overrideEnergy; }

  G4int remainingIndex = electronIndex;
  ShellSpec shell = kShells.back();
  for(const auto& candidate : kShells) {
    if(remainingIndex < candidate.capacity) { shell = candidate; break; }
    remainingIndex -= candidate.capacity;
  }

  const G4double innerShield = 0.30*static_cast<G4double>(std::max(0, electronIndex));
  const G4double zEff = std::max(1.0, static_cast<G4double>(Zp) - innerShield);
  const G4double rydberg = 13.605693122994*eV;
  return rydberg*zEff*zEff/(static_cast<G4double>(shell.n*shell.n))*shell.relFactor;
}

G4CCBindingEnergyProviderStore& G4CCBindingEnergyProviderStore::Instance()
{
  static G4CCBindingEnergyProviderStore store;
  return store;
}

G4CCBindingEnergyProviderStore::G4CCBindingEnergyProviderStore()
  : fProvider(std::make_unique<G4CCHydrogenicBindingEnergyProvider>())
{}

const G4CCBindingEnergyProvider& G4CCBindingEnergyProviderStore::Provider() const
{
  return *fProvider;
}

void G4CCBindingEnergyProviderStore::SetProvider(std::unique_ptr<G4CCBindingEnergyProvider> provider)
{
  if(provider) { fProvider = std::move(provider); }
}
