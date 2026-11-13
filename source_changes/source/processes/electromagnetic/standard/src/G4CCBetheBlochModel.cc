#include "G4CCBetheBlochModel.hh"

#include "G4CCIonChargeState.hh"
#include "G4CCCurrentChargeConfig.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicalConstants.hh"
#include "G4IonisParamMat.hh"
#include "G4NistManager.hh"
#include "G4LossTableManager.hh"
#include "G4Log.hh"
#include "G4EmCorrections.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VEmFluctuationModel.hh"

#include <algorithm>
#include <cmath>

namespace
{
G4bool G4CCUseLowEnergyVanillaFallback(const G4ParticleDefinition* particle,
                                       const G4Material* material,
                                       G4double kineticEnergy)
{
  if(!G4CCIonChargeState::IonisationModelApplies(particle, material)) { return false; }
  const G4double eMeVPerU =
    G4CCIonChargeState::KineticEnergyPerNucleon(particle, kineticEnergy)/CLHEP::MeV;
  return eMeVPerU <= (G4CCCurrentChargeConfig::Instance().IsDiracMode()
    ? G4CCCurrentChargeConfig::Low0MeVPerU() : 30.0);
}
}

G4CCBetheBlochModel::G4CCBetheBlochModel(const G4ParticleDefinition* p,
                                         const G4String& name)
  : G4BetheBlochModel(p, name)
{}

G4double G4CCBetheBlochModel::ChargeSquareRatio(const G4Track& track)
{
  if(!G4CCIonChargeState::IonisationModelApplies(track) ||
     G4CCUseLowEnergyVanillaFallback(track.GetParticleDefinition(),
                                     track.GetMaterial(),
                                     track.GetKineticEnergy())) {
    fHasCachedCurrentCharge = false;
    return G4BetheBlochModel::ChargeSquareRatio(track);
  }

  fCachedChargeInEplus = G4CCIonChargeState::GetChargeInEplus(track);
  fHasCachedCurrentCharge = true;

  const G4double q2 = std::max(0.0, G4CCIonChargeState::StoppingChargeSquare(track));
  if(auto* fluct = GetModelOfFluctuations()) {
    fluct->SetParticleAndCharge(track.GetParticleDefinition(), q2);
  }
  return q2;
}

G4double G4CCBetheBlochModel::CrossSectionPerVolume(const G4Material* material,
                                                    const G4ParticleDefinition* particle,
                                                    G4double kineticEnergy,
                                                    G4double cutEnergy,
                                                    G4double maxEnergy)
{
  if(G4CCIonChargeState::IonisationModelApplies(particle, material) &&
     G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(
       G4CCIonChargeState::KineticEnergyPerNucleon(particle, kineticEnergy)/CLHEP::MeV)) {
    return 0.0;
  }
  return G4BetheBlochModel::CrossSectionPerVolume(material, particle, kineticEnergy,
                                                  cutEnergy, maxEnergy);
}

G4double G4CCBetheBlochModel::ComputeDEDXPerVolume(const G4Material* material,
                                                   const G4ParticleDefinition* particle,
                                                   G4double kineticEnergy,
                                                   G4double cut)
{
  if(!G4CCIonChargeState::IonisationModelApplies(particle, material) ||
     G4CCUseLowEnergyVanillaFallback(particle, material, kineticEnergy)) {
    return G4BetheBlochModel::ComputeDEDXPerVolume(material, particle, kineticEnergy, cut);
  }

  const G4double mass = particle->GetPDGMass();
  if(mass <= 0.0 || kineticEnergy <= 0.0 || material == nullptr) { return 0.0; }

  const G4double tmax = MaxSecondaryEnergy(particle, kineticEnergy);

  const G4double spin = particle->GetPDGSpin();
  const G4double qDefinition = particle->GetPDGCharge()*inveplus;
  G4int zForFormFactor = G4CCIonChargeState::GetZp(particle);
  if(zForFormFactor <= 0) { zForFormFactor = G4lrint(std::abs(qDefinition)); }

  auto* nist = G4NistManager::Instance();
  G4double xff = 0.8426*CLHEP::GeV;
  if(spin == 0.0 && mass < CLHEP::GeV) {
    xff = 0.736*CLHEP::GeV;
  } else if(mass > CLHEP::GeV && zForFormFactor > 1) {
    xff /= nist->GetA27(zForFormFactor);
  }
  const G4double formfactor = 2.0*CLHEP::electron_mass_c2/(xff*xff);
  const G4double tlimit = 2.0/formfactor;
  const G4double cutEnergy = std::min(std::min(cut, tmax), tlimit);
  if(cutEnergy <= 0.0 || tmax <= 0.0) { return 0.0; }

  const G4double tau = kineticEnergy/mass;
  const G4double gamma = tau + 1.0;
  const G4double bg2 = tau*(tau + 2.0);
  const G4double beta2 = bg2/(gamma*gamma);
  if(beta2 <= 0.0) { return 0.0; }

  const G4double xc = cutEnergy/tmax;
  const G4double meanExcitation = material->GetIonisation()->GetMeanExcitationEnergy();
  const G4double meanExcitation2 = meanExcitation*meanExcitation;

  auto* corr = G4LossTableManager::Instance()->EmCorrections();

  const G4double q2Table = GetChargeSquareRatio(particle, material, kineticEnergy);

  G4double dedx = G4Log(2.0*CLHEP::electron_mass_c2*bg2*cutEnergy/meanExcitation2)
                - (1.0 + xc)*beta2;

  if(spin > 0.0) {
    const G4double del = 0.5*cutEnergy/(kineticEnergy + mass);
    dedx += del*del;
  }

  const G4double twoln10 = 2.0*G4Log(10.0);
  const G4double x = G4Log(bg2)/twoln10;
  dedx -= material->GetIonisation()->DensityCorrection(x);
  dedx -= 2.0*corr->ShellCorrection(particle, material, kineticEnergy);

  dedx *= CLHEP::twopi_mc2_rcl2*q2Table*material->GetElectronDensity()/beta2;
  return std::max(dedx, 0.0);
}

G4double G4CCBetheBlochModel::GetChargeSquareRatio(const G4ParticleDefinition* particle,
                                                   const G4Material* material,
                                                   G4double kineticEnergy)
{
  if(G4CCIonChargeState::IonisationModelApplies(particle, material) &&
     !G4CCUseLowEnergyVanillaFallback(particle, material, kineticEnergy)) {
    return 1.0;
  }
  return G4BetheBlochModel::GetChargeSquareRatio(particle, material, kineticEnergy);
}

G4double G4CCBetheBlochModel::GetParticleCharge(const G4ParticleDefinition* particle,
                                                const G4Material* material,
                                                G4double kineticEnergy)
{
  if(G4CCIonChargeState::IonisationModelApplies(particle, material) &&
     !G4CCUseLowEnergyVanillaFallback(particle, material, kineticEnergy)) {
    if(fHasCachedCurrentCharge) { return fCachedChargeInEplus*eplus; }
    return static_cast<G4double>(G4CCIonChargeState::GetZp(particle))*eplus;
  }
  return G4BetheBlochModel::GetParticleCharge(particle, material, kineticEnergy);
}

void G4CCBetheBlochModel::CorrectionsAlongStep(const G4Material* material,
                                               const G4ParticleDefinition* particle,
                                               const G4double kineticEnergy,
                                               const G4double cutEnergy,
                                               const G4double& length,
                                               G4double& eloss)
{
  if(G4CCIonChargeState::IonisationModelApplies(particle, material) &&
     !G4CCUseLowEnergyVanillaFallback(particle, material, kineticEnergy)) {
    return;
  }
  G4BetheBlochModel::CorrectionsAlongStep(material, particle, kineticEnergy,
                                          cutEnergy, length, eloss);
}

void G4CCBetheBlochModel::SampleSecondaries(std::vector<G4DynamicParticle*>* secondaries,
                                            const G4MaterialCutsCouple* couple,
                                            const G4DynamicParticle* dynamicParticle,
                                            G4double tmin,
                                            G4double maxEnergy)
{
  const auto* material = couple == nullptr ? nullptr : couple->GetMaterial();
  if(dynamicParticle != nullptr &&
     G4CCIonChargeState::IonisationModelApplies(dynamicParticle->GetDefinition(), material) &&
     G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(
       G4CCIonChargeState::KineticEnergyPerNucleon(dynamicParticle->GetDefinition(),
                                                 dynamicParticle->GetKineticEnergy())/CLHEP::MeV)) {
    return;
  }
  G4BetheBlochModel::SampleSecondaries(secondaries, couple, dynamicParticle, tmin, maxEnergy);
}
