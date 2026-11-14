#include "G4CCScreening.hh"

#include "G4CCBindingEnergyProvider.hh"
#include "G4CCIonChargeState.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"

#include <algorithm>
#include <cmath>

G4double G4CCScreening::fExponent = 2.0;

G4double G4CCScreening::Exponent()
{
  return fExponent;
}

void G4CCScreening::SetExponent(G4double exponent)
{
  if(exponent > 0.0) { fExponent = exponent; }
}

G4double G4CCScreening::ScreeningFunction(G4double x)
{
  if(x <= 0.0) { return 1.0; }
  return 1.0/(1.0 + std::pow(x, fExponent));
}

G4double G4CCScreening::DeltaCharge(G4int Zp, G4int Qe, G4double deltaElectronEnergy)
{
  if(Zp <= 0) { return 0.0; }
  Qe = G4CCIonChargeState::ClampQe(Zp, Qe);
  if(Qe == 0) { return static_cast<G4double>(Zp); }

  const auto& provider = G4CCBindingEnergyProviderStore::Instance().Provider();
  G4double screening = 0.0;
  for(G4int i = 0; i < Qe; ++i) {
    const G4double Ii = provider.BindingEnergy(Zp, i);
    if(Ii <= 0.0) { continue; }
    screening += ScreeningFunction(deltaElectronEnergy/Ii);
  }
  const G4double q = static_cast<G4double>(Zp) - screening;
  return std::max(0.0, q);
}

G4double G4CCScreening::DeltaChargeSquare(G4int Zp, G4int Qe, G4double deltaElectronEnergy)
{
  const G4double q = DeltaCharge(Zp, Qe, deltaElectronEnergy);
  return q*q;
}

G4double G4CCScreening::DeltaChargeSquare(const G4DynamicParticle* projectile,
                                          G4double deltaElectronEnergy)
{
  if(projectile == nullptr) { return 0.0; }
  const auto* particle = projectile->GetDefinition();
  const G4int Zp = G4CCIonChargeState::GetZp(particle);
  const G4int Qe = G4CCIonChargeState::GetQe(particle,
                                             G4CCIonChargeState::GetChargeInEplus(projectile));
  return DeltaChargeSquare(Zp, Qe, deltaElectronEnergy);
}
