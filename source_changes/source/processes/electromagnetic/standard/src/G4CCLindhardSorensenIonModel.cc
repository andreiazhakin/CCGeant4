#include "G4CCLindhardSorensenIonModel.hh"

#include "G4AutoLock.hh"
#include "G4BetheBlochModel.hh"
#include "G4CCAtimaQDedxTable.hh"
#include "G4CCAtimaEquilibriumDedxTable.hh"
#include "G4CCBetheBlochModel.hh"
#include "G4CCBraggIonModel.hh"
#include "G4CCCurrentChargeConfig.hh"
#include "G4CCIonChargeState.hh"
#include "G4CCMottCorrection.hh"
#include "G4CCIonDeltaElectronProcess.hh"
#include "G4CCScreening.hh"
#include "G4CCDiracPhysics.hh"
#include "G4CCDiracTable.hh"
#include "G4DeltaAngle.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4EmCorrections.hh"
#include "G4Exception.hh"
#include "G4IonICRU73Data.hh"
#include "G4LindhardSorensenData.hh"
#include "G4Log.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VEmFluctuationModel.hh"
#include "Randomize.hh"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <string>

G4LindhardSorensenData* G4CCLindhardSorensenIonModel::fLSData = nullptr;
G4IonICRU73Data* G4CCLindhardSorensenIonModel::fIonData = nullptr;


namespace
{
thread_local const G4Track* gCCLSTransportTrack = nullptr;

G4bool G4CCLSEnvFlagValue(const char* value, G4bool defaultValue)
{
  if(value == nullptr) { return defaultValue; }
  const std::string s(value);
  if(s == "0" || s == "false" || s == "FALSE" ||
     s == "off" || s == "OFF" || s == "no" || s == "NO") {
    return false;
  }
  if(s == "1" || s == "true" || s == "TRUE" ||
     s == "on" || s == "ON" || s == "yes" || s == "YES") {
    return true;
  }
  return defaultValue;
}

G4bool G4CCLSEnvFlag(const char* name, G4bool defaultValue)
{
  return G4CCLSEnvFlagValue(std::getenv(name), defaultValue);
}

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

G4bool G4CCLSUseMottDeltaRuntime()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return true;
  static const char* envValue = std::getenv("G4CC_USE_MOTT_DELTA");
  static const G4bool envIsSet = (envValue != nullptr);
  static const G4bool cachedEnv = G4CCLSEnvFlagValue(envValue, true);
  if(envIsSet) { return cachedEnv; }
  return G4CCCurrentChargeConfig::Instance().UseMottDeltaElectrons();
}

G4bool G4CCLSDisableContinuousLSRuntime()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return true;
  static const char* envValue = std::getenv("G4CC_DISABLE_CONTINUOUS_LS_FOR_MOTT_DELTA");
  static const G4bool envIsSet = (envValue != nullptr);
  static const G4bool cachedEnv = G4CCLSEnvFlagValue(envValue, true);
  if(envIsSet) { return cachedEnv; }
  return G4CCCurrentChargeConfig::Instance().DisableContinuousLSWhenMottDelta();
}

G4bool G4CCLSDiagnosticsEnabled()
{
  static const G4bool enabled = G4CCLSEnvFlag("G4CC_MOTT_DIAGNOSTICS", false);
  return enabled;
}

G4bool G4CCDirectBlochEnabled()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return true;
  static const G4bool enabled = G4CCLSEnvFlag("G4CC_ENABLE_DIRECT_BLOCH", true);
  return enabled;
}

G4bool G4CCBarkasEnabled()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return true;
  static const G4bool enabled = G4CCLSEnvFlag("G4CC_ENABLE_BARKAS_L1", true);
  return enabled;
}

G4bool G4CCGeant4HighOrderFEnabled()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return false;
  static const G4bool enabled = G4CCLSEnvFlag("G4CC_ENABLE_GEANT4_HIGH_ORDER_F", false);
  return enabled;
}

G4bool G4CCBlochDiagnosticsEnabled()
{
  static const G4bool enabled = G4CCLSEnvFlag("G4CC_BLOCH_DIAGNOSTICS", false);
  return enabled;
}

G4int G4CCDirectBlochTerms()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return 1000;
  static const G4int terms = []() {
    const char* value = std::getenv("G4CC_DIRECT_BLOCH_TERMS");
    if(value == nullptr) { return 1000; }
    const G4int n = std::atoi(value);
    return std::max(1, n);
  }();
  return terms;
}

G4String G4CCDefaultCorrectionZMode()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return G4String("effective");
  static const G4String mode = []() {
    const char* value = std::getenv("G4CC_CORRECTION_Z_MODE");
    if(value == nullptr) { value = std::getenv("G4CC_CORRECTION_CHARGE_MODE"); }
    if(value == nullptr) { return G4String("effective"); }
    return G4String(value);
  }();
  return mode;
}

G4String G4CCBlochZMode()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return G4String("effective");
  static const G4String mode = []() {
    const char* value = std::getenv("G4CC_BLOCH_Z_MODE");
    if(value == nullptr) { return G4CCDefaultCorrectionZMode(); }
    return G4String(value);
  }();
  return mode;
}

G4String G4CCBarkasZMode()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return G4String("effective");
  static const G4String mode = []() {
    const char* value = std::getenv("G4CC_BARKAS_Z_MODE");
    if(value == nullptr) { return G4CCDefaultCorrectionZMode(); }
    return G4String(value);
  }();
  return mode;
}

G4double G4CCDirectBlochDeltaL(G4double zCorr, G4double beta)
{
  if(zCorr == 0.0 || beta <= 0.0) { return 0.0; }

  const G4double alpha = CLHEP::fine_structure_const;
  const G4double eta = std::abs(zCorr)*alpha/beta;
  const G4double eta2 = eta*eta;
  const G4int nTerms = G4CCDirectBlochTerms();

  G4double sum = 0.0;
  for(G4int n = 1; n <= nTerms; ++n) {
    const G4double x = static_cast<G4double>(n);
    sum += 1.0/(x*(x*x + eta2));
  }

  return -eta2*sum;
}

G4double G4CCChooseCorrectionCharge(const G4String& mode,
                                    G4double currentChargeInEplus,
                                    G4int zProjectile,
                                    G4EmCorrections* corrections,
                                    const G4ParticleDefinition* p,
                                    const G4Material* mat,
                                    G4double kineticEnergy)
{
  if(mode == "projectile" || mode == "Zp" || mode == "zp" || mode == "nuclear") {
    return static_cast<G4double>(std::max(1, zProjectile));
  }

  if(mode == "effective" || mode == "eff" || mode == "qeff") {
    if(corrections != nullptr && p != nullptr && mat != nullptr) {
      return std::abs(corrections->GetParticleCharge(p, mat, kineticEnergy)/CLHEP::eplus);
    }
  }

  if(currentChargeInEplus > 0.0) { return currentChargeInEplus; }
  return static_cast<G4double>(std::max(1, zProjectile));
}


std::atomic<int> gLSSuppressionDiagnosticPrints{0};
std::atomic<int> gBlochDiagnosticPrints{0};

std::atomic<int> gAtimaDeltaSubDiagnosticPrints{0};

G4bool G4CCAtimaDeltaSubDiagnosticsEnabled()
{
  static const G4bool enabled = G4CCLSEnvFlag("G4CC_ATIMA_DELTA_SUB_DIAGNOSTICS", false);
  return enabled;
}

G4String G4CCMottZModeForDeltaSubtraction()
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return G4String("effective");
  static const G4String mode = []() {
    const char* value = std::getenv("G4CC_MOTT_Z_MODE");
    if(value == nullptr) { return G4CCDefaultCorrectionZMode(); }
    return G4String(value);
  }();
  return mode;
}

G4int G4CCChooseMottZForDeltaSubtraction(G4double currentChargeInEplus,
                                          G4int zProjectile,
                                          G4EmCorrections* corrections,
                                          const G4ParticleDefinition* p,
                                          const G4Material* mat,
                                          G4double kineticEnergy)
{
  if(zProjectile <= 0) { return 0; }

  const G4double z = G4CCChooseCorrectionCharge(G4CCMottZModeForDeltaSubtraction(),
                                                currentChargeInEplus, zProjectile,
                                                corrections, p, mat, kineticEnergy);
  G4int zInt = std::max(1, static_cast<G4int>(std::lround(z)));
  if(!G4CCMottCorrection::IsAvailable(zInt) &&
     G4CCMottCorrection::IsAvailable(zProjectile)) {
    zInt = zProjectile;
  }
  return zInt;
}

G4double G4CCFiniteSizeLossIntegral(const G4ParticleDefinition* particle,
                                   G4double kineticEnergy, G4int zMott)
{
  if(particle == nullptr || kineticEnergy <= 0.0 ||
     !G4CCMottCorrection::IsAvailable(zMott)) { return 0.0; }
  const G4double threshold =
    G4CCIonDeltaElectronProcess::FormFactorThresholdEnergy(particle, kineticEnergy);
  const G4double mass = particle->GetPDGMass();
  if(mass <= 0.0) { return 0.0; }
  const G4double tau = kineticEnergy/mass;
  const G4double gamma = 1.0 + tau;
  const G4double ratio = CLHEP::electron_mass_c2/mass;
  const G4double tMax = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.0)
    /(1.0 + 2.0*gamma*ratio + ratio*ratio);
  if(!(threshold > 0.0 && threshold < tMax)) { return 0.0; }

  const G4double spin = particle->GetPDGSpin();
  const G4double magneticMoment = particle->GetPDGMagneticMoment();
  struct Cache {
    const G4ParticleDefinition* particle = nullptr;
    G4double energy = 0.0, mass = 0.0, spin = 0.0, magneticMoment = 0.0;
    G4double threshold = 0.0, integral = 0.0;
    G4int zMott = 0;
  };
  static thread_local Cache cache;
  if(cache.particle == particle && cache.energy == kineticEnergy &&
     cache.mass == mass && cache.spin == spin &&
     cache.magneticMoment == magneticMoment && cache.threshold == threshold &&
     cache.zMott == zMott) { return cache.integral; }

  constexpr G4double nodes[4] = {0.1834346424956498, 0.5255324099163290,
                                0.7966664774136267, 0.9602898564975363};
  constexpr G4double weights[4] = {0.3626837833783620, 0.3137066458778873,
                                  0.2223810344533745, 0.1012285362903763};
  const G4double beta = std::sqrt(tau*(tau + 2.0)/(gamma*gamma));
  const G4double low = std::log(threshold);
  const G4double span = std::log(tMax/threshold);
  const G4int panels = std::max(1, static_cast<G4int>(std::ceil(2.0*span)));
  const G4double du = span/panels;
  G4double integral = 0.0;
  for(G4int panel = 0; panel < panels; ++panel) {
    const G4double mid = low + (panel + 0.5)*du;
    for(G4int j = 0; j < 4; ++j) for(G4int sign : {-1, 1}) {
      const G4double t = std::exp(mid + sign*0.5*du*nodes[j]);
      const G4double rejected = 1.0 -
        G4CCIonDeltaElectronProcess::FormFactorAcceptance(particle, kineticEnergy, t);
      integral += 0.5*du*weights[j]*rejected*
        G4CCMottCorrection::RMott(zMott, beta, t, tMax);
    }
  }
  cache = {particle, kineticEnergy, mass, spin, magneticMoment,
           threshold, std::max(0.0, integral), zMott};
  return cache.integral;
}

G4double G4CCDeltaMeanLossIntegralLog(G4int zProjectile,
                                      G4int zMott,
                                      G4int qe,
                                      G4double tMin,
                                      G4double tMax,
                                      G4double beta2,
                                      G4double totalEnergy2,
                                      G4double spin,
                                      G4double qEffective2,
                                      G4double globalWeight,
                                      G4bool useMott,
                                      const G4ParticleDefinition* particle,
                                      G4double kineticEnergy,
                                      G4bool applyFallbackFormFactor)
{
  if(tMax <= tMin || tMin <= 0.0 || beta2 <= 0.0 || zProjectile <= 0) { return 0.0; }

  const G4double beta = std::sqrt(beta2);
  if(applyFallbackFormFactor && particle != nullptr) {
    constexpr G4double nodes[4] = {0.1834346424956498, 0.5255324099163290,
                                  0.7966664774136267, 0.9602898564975363};
    constexpr G4double weights[4] = {0.3626837833783620, 0.3137066458778873,
                                    0.2223810344533745, 0.1012285362903763};
    const G4double split = G4CCIonDeltaElectronProcess::FormFactorThresholdEnergy(particle, kineticEnergy);
    const G4double boundary = std::max(tMin, std::min(tMax, split));
    const G4double bounds[3] = {tMin, boundary, tMax};
    const G4double w = std::max(0.0, std::min(1.0, globalWeight));
    G4double acceptedIntegral = 0.0;
    for(G4int segment = 0; segment < 2; ++segment) {
      if(bounds[segment+1] <= bounds[segment]) { continue; }
      const G4double low = std::log(bounds[segment]);
      const G4double du = std::log(bounds[segment+1]/bounds[segment])/48.0;
      for(G4int panel = 0; panel < 48; ++panel) {
        const G4double mid = low + (panel + 0.5)*du;
        for(G4int j = 0; j < 4; ++j) for(G4int sign : {-1, 1}) {
          const G4double t = std::exp(mid + sign*0.5*du*nodes[j]);
          const G4double current = w > 0.0 ? G4CCScreening::DeltaChargeSquare(zProjectile, qe, t) : 0.0;
          const G4double ff = G4CCIonDeltaElectronProcess::FormFactorAcceptance(particle, kineticEnergy, t);
          const G4double charge2 = w*current + (1.0-w)*std::max(0.0, qEffective2)*ff;
          const G4double shape = useMott
            ? G4CCMottCorrection::RMott(zMott, beta, t, tMax)
            : G4CCMottCorrection::G4BornShape(beta2, t, tMax, totalEnergy2, spin);
          acceptedIntegral += 0.5*du*weights[j]*charge2*shape;
        }
      }
    }
    return std::max(0.0, acceptedIntegral);
  }
  constexpr G4int nBins = 96;
  const G4double u0 = std::log(tMin);
  const G4double u1 = std::log(tMax);
  const G4double du = (u1 - u0)/static_cast<G4double>(nBins);

  G4double integral = 0.0;
  for(G4int i = 0; i < nBins; ++i) {
    const G4double ua = u0 + du*static_cast<G4double>(i);
    const G4double ub = ua + du;
    const G4double um = 0.5*(ua + ub);
    const G4double ta = std::exp(ua);
    const G4double tm = std::exp(um);
    const G4double tb = std::exp(ub);
    const auto f = [&](G4double t) {
      const G4double w = std::max(0.0, std::min(1.0, globalWeight));
      const G4double z2Delta = w > 0.0
        ? G4CCScreening::DeltaChargeSquare(zProjectile, qe, t) : 0.0;
      const G4double z2 = std::max(0.0, w*z2Delta + (1.0 - w)*std::max(0.0, qEffective2));
      const G4double shape = useMott
        ? G4CCMottCorrection::RMott(zMott, beta, t, tMax)
        : G4CCMottCorrection::G4BornShape(beta2, t, tMax, totalEnergy2, spin);
      return z2*shape;
    };
    integral += (du/6.0)*(f(ta) + 4.0*f(tm) + f(tb));
  }
  return std::max(0.0, integral);
}

G4double G4CCDeltaMeanLossPerVolumeForChargeState(const G4Material* mat,
                                                  const G4ParticleDefinition* p,
                                                  G4double kineticEnergy,
                                                  G4double cutEnergy,
                                                  G4int zProjectile,
                                                  G4int qe,
                                                  G4double currentChargeInEplus,
                                                  G4double qEffective2,
                                                  G4double globalWeight,
                                                  G4bool mottActive,
                                                  G4EmCorrections* corrections,
                                                  G4bool applyFallbackFormFactor = false)
{
  if(mat == nullptr || p == nullptr || kineticEnergy <= 0.0 || cutEnergy <= 0.0) { return 0.0; }

  const G4double mass = p->GetPDGMass();
  if(mass <= 0.0) { return 0.0; }

  const G4double tau = kineticEnergy/mass;
  const G4double electronMassRatio = CLHEP::electron_mass_c2/mass;
  const G4double tMax = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.0)
    /(1.0 + 2.0*(tau + 1.0)*electronMassRatio + electronMassRatio*electronMassRatio);
  const G4double tMin = std::min(cutEnergy, tMax);
  if(tMax <= tMin || tMin <= 0.0) { return 0.0; }

  const G4double totalEnergy = kineticEnergy + mass;
  const G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*mass)/(totalEnergy*totalEnergy);
  if(beta2 <= 0.0) { return 0.0; }

  const G4int zMott = G4CCChooseMottZForDeltaSubtraction(currentChargeInEplus,
                                                         zProjectile, corrections,
                                                         p, mat, kineticEnergy);
  const G4bool useMott = mottActive && G4CCMottCorrection::IsAvailable(zMott);

  const G4double integral = G4CCDeltaMeanLossIntegralLog(zProjectile, zMott, qe,
                                                        tMin, tMax, beta2,
                                                        totalEnergy*totalEnergy,
                                                        p->GetPDGSpin(),
                                                        qEffective2, globalWeight,
                                                        useMott, p, kineticEnergy,
                                                        applyFallbackFormFactor);
  if(integral <= 0.0) { return 0.0; }
  return mat->GetElectronDensity()*CLHEP::twopi_mc2_rcl2*integral/beta2;
}

}

namespace
{
  G4Mutex ccIonXSMutex = G4MUTEX_INITIALIZER;
}

G4CCLindhardSorensenIonModel::G4CCLindhardSorensenIonModel(
  const G4ParticleDefinition*, const G4String& name)
  : G4VEmModel(name), fTwoLog10(2.0*G4Log(10.0))
{
  fElectron = G4Electron::Electron();
  fCorrections = G4LossTableManager::Instance()->EmCorrections();
  fNist = G4NistManager::Instance();

  fBraggModel = new G4CCBraggIonModel();
  fBetheBlochModel = new G4CCBetheBlochModel();
  fLowEnergyLimit = 2.0*CLHEP::MeV;
}

G4CCLindhardSorensenIonModel::~G4CCLindhardSorensenIonModel()
{
  delete fBraggModel;
  delete fBetheBlochModel;

  if(fOwnsSharedData) {
    delete fLSData;
    delete fIonData;
    fLSData = nullptr;
    fIonData = nullptr;
  }
}

void G4CCLindhardSorensenIonModel::Initialise(const G4ParticleDefinition* p,
                                              const G4DataVector& cuts)
{
  fBraggModel->Initialise(p, cuts);
  fBetheBlochModel->Initialise(p, cuts);
  SetParticle(p);

  SetDeexcitationFlag(false);

  if(fParticleChangeForLoss == nullptr) {
    fParticleChangeForLoss = GetParticleChangeForLoss();
    if(UseAngularGeneratorFlag() && GetAngularDistribution() == nullptr) {
      SetAngularDistribution(new G4DeltaAngle());
    }
  }

  if(fLSData == nullptr) {
    G4AutoLock lock(&ccIonXSMutex);
    if(fLSData == nullptr) {
      fOwnsSharedData = true;
      fLSData = new G4LindhardSorensenData();
      fIonData = new G4IonICRU73Data();
    }
    lock.unlock();
  }

  if(fOwnsSharedData) {
    fIonData->Initialise();
  }
}

G4CCLindhardSorensenIonModel::TransportScope::TransportScope(const G4Track& track)
  : fPrevious(gCCLSTransportTrack)
{
  gCCLSTransportTrack = &track;
}

G4CCLindhardSorensenIonModel::TransportScope::~TransportScope()
{
  gCCLSTransportTrack = fPrevious;
}

G4bool G4CCLindhardSorensenIonModel::InTransportContext(
  const G4ParticleDefinition* p, const G4Material* mat)
{
  return gCCLSTransportTrack != nullptr &&
    gCCLSTransportTrack->GetParticleDefinition() == p &&
    gCCLSTransportTrack->GetMaterial() == mat;
}

G4bool G4CCLindhardSorensenIonModel::MayUseCachedCharge(
  const G4ParticleDefinition* p, const G4Material* mat) const
{
  return fHasCachedCurrentCharge && InTransportContext(p, mat);
}

G4double G4CCLindhardSorensenIonModel::ChargeSquareRatio(const G4Track& track)
{
  if(G4CCCurrentChargeConfig::Instance().UseAtimaEquilibriumDedx() &&
     G4CCIonChargeState::IonisationModelApplies(track)) {
    fCachedChargeInEplus = G4CCIonChargeState::GetChargeInEplus(track);
    fHasCachedCurrentCharge = true;
  }
  if(!G4CCIonChargeState::IonisationModelApplies(track) ||
     G4CCUseLowEnergyVanillaFallback(track.GetParticleDefinition(),
                                     track.GetMaterial(),
                                     track.GetKineticEnergy())) {
    if(!G4CCCurrentChargeConfig::Instance().UseAtimaEquilibriumDedx() ||
       !G4CCIonChargeState::IonisationModelApplies(track)) {
      fHasCachedCurrentCharge = false;
    }
    return GetChargeSquareRatio(track.GetParticleDefinition(),
                                track.GetMaterial(),
                                track.GetKineticEnergy());
  }

  fCachedChargeInEplus = G4CCIonChargeState::GetChargeInEplus(track);
  fHasCachedCurrentCharge = true;

  fChargeSquareRatio = std::max(0.0, G4CCIonChargeState::StoppingChargeSquare(track));
  fBraggModel->SetChargeSquareRatio(fChargeSquareRatio);
  fBetheBlochModel->SetChargeSquareRatio(fChargeSquareRatio);

  if(auto* fluct = GetModelOfFluctuations()) {
    fluct->SetParticleAndCharge(track.GetParticleDefinition(), fChargeSquareRatio);
  }
  return fChargeSquareRatio;
}

G4double G4CCLindhardSorensenIonModel::GetChargeSquareRatio(
  const G4ParticleDefinition* p, const G4Material* mat, G4double kineticEnergy)
{
  if(G4CCIonChargeState::IonisationModelApplies(p, mat) &&
     !G4CCUseLowEnergyVanillaFallback(p, mat, kineticEnergy)) {
    fChargeSquareRatio = 1.0;
  } else {
    fCorrections->GetParticleCharge(p, mat, kineticEnergy);
    fChargeSquareRatio = fCorrections->EffectiveChargeSquareRatio(p, mat, kineticEnergy);
  }

  fBraggModel->SetChargeSquareRatio(fChargeSquareRatio);
  fBetheBlochModel->SetChargeSquareRatio(fChargeSquareRatio);
  return fChargeSquareRatio;
}

G4double G4CCLindhardSorensenIonModel::GetParticleCharge(
  const G4ParticleDefinition* p, const G4Material* mat, G4double kineticEnergy)
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode() &&
     G4CCIonChargeState::IonisationModelApplies(p, mat)) {
    const G4double eu = G4CCIonChargeState::KineticEnergyPerNucleon(p, kineticEnergy)/MeV;
    if(eu <= G4CCCurrentChargeConfig::Low0MeVPerU()) {
      return fCorrections->GetParticleCharge(p, mat, kineticEnergy);
    }
    if(eu >= G4CCCurrentChargeConfig::High1MeVPerU()) {
      return static_cast<G4double>(G4CCIonChargeState::GetZp(p))*eplus;
    }
    if(MayUseCachedCharge(p, mat)) { return fCachedChargeInEplus*eplus; }
    return static_cast<G4double>(G4CCIonChargeState::GetZp(p))*eplus;
  }
  if(G4CCIonChargeState::IonisationModelApplies(p, mat) &&
     !G4CCUseLowEnergyVanillaFallback(p, mat, kineticEnergy)) {
    if(MayUseCachedCharge(p, mat)) { return fCachedChargeInEplus*eplus; }
    return static_cast<G4double>(G4CCIonChargeState::GetZp(p))*eplus;
  }
  return fCorrections->GetParticleCharge(p, mat, kineticEnergy);
}

void G4CCLindhardSorensenIonModel::SetupParameters()
{
  fMass = fParticle->GetPDGMass();
  fSpin = fParticle->GetPDGSpin();
  fDefinitionChargeInEplus = fParticle->GetPDGCharge()*inveplus;

  fZin = G4lrint(std::abs(fDefinitionChargeInEplus));
  const G4int zProjectile = G4CCIonChargeState::GetZp(fParticle);
  if(zProjectile > 0) { fZin = zProjectile; }
  fZin = std::max(1, fZin);

  fElectronMassRatio = CLHEP::electron_mass_c2/fMass;
  fProtonMassRatio = CLHEP::proton_mass_c2/fMass;

  const G4double aMag = 1.0/(0.5*CLHEP::eplus*CLHEP::hbar_Planck*CLHEP::c_squared);
  const G4double magmom = fParticle->GetPDGMagneticMoment()*fMass*aMag;
  fMagMoment2 = magmom*magmom - 1.0;

  G4double x = 0.8426*CLHEP::GeV;
  if(fSpin == 0.0 && fMass < CLHEP::GeV) {
    x = 0.736*CLHEP::GeV;
  } else if(fZin > 1) {
    x /= fNist->GetA27(fZin);
  }

  fFormFactor = 2.0*CLHEP::electron_mass_c2/(x*x);
  fTlimit = 2.0/fFormFactor;
}

G4double G4CCLindhardSorensenIonModel::MinEnergyCut(
  const G4ParticleDefinition*, const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

G4double G4CCLindhardSorensenIonModel::ComputeCrossSectionPerElectron(
  const G4ParticleDefinition* p, G4double kineticEnergy,
  G4double cutEnergy, G4double maxKinEnergy)
{
  const G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  const G4double emax = std::min(tmax, maxKinEnergy);
  const G4double escaled = kineticEnergy*fProtonMassRatio;

  return (escaled <= fLowEnergyLimit)
    ? fBraggModel->ComputeCrossSectionPerElectron(p, kineticEnergy, cutEnergy, emax)
    : fBetheBlochModel->ComputeCrossSectionPerElectron(p, kineticEnergy, cutEnergy, emax);
}

G4double G4CCLindhardSorensenIonModel::ComputeCrossSectionPerAtom(
  const G4ParticleDefinition* p, G4double kineticEnergy, G4double Z, G4double,
  G4double cutEnergy, G4double maxEnergy)
{
  return Z*ComputeCrossSectionPerElectron(p, kineticEnergy, cutEnergy, maxEnergy);
}

G4double G4CCLindhardSorensenIonModel::CrossSectionPerVolume(
  const G4Material* material, const G4ParticleDefinition* p,
  G4double kineticEnergy, G4double cutEnergy, G4double maxEnergy)
{
  if(G4CCIonChargeState::IonisationModelApplies(p, material) &&
     G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(
       G4CCIonChargeState::KineticEnergyPerNucleon(p, kineticEnergy)/CLHEP::MeV)) {
    return 0.0;
  }

  GetChargeSquareRatio(p, material, kineticEnergy);
  return material->GetElectronDensity()
    *ComputeCrossSectionPerElectron(p, kineticEnergy, cutEnergy, maxEnergy);
}

G4double G4CCLindhardSorensenIonModel::ComputeOldUnitDEDXPerVolume(
  const G4Material* mat, const G4ParticleDefinition* p,
  G4double kineticEnergy, G4double cut)
{
  GetChargeSquareRatio(p, mat, kineticEnergy);

  const G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  const G4double cutEnergy = std::min(std::min(cut, tmax), fTlimit);
  const G4double escaled = kineticEnergy*fProtonMassRatio;

  return (escaled <= fLowEnergyLimit)
    ? fBraggModel->ComputeDEDXPerVolume(mat, p, kineticEnergy, cutEnergy)
    : fBetheBlochModel->ComputeDEDXPerVolume(mat, p, kineticEnergy, cutEnergy);
}

G4double G4CCLindhardSorensenIonModel::ComputeDEDXPerVolume(
  const G4Material* mat, const G4ParticleDefinition* p,
  G4double kineticEnergy, G4double cut)
{
  const G4double oldUnitDedx = ComputeOldUnitDEDXPerVolume(mat, p, kineticEnergy, cut);
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode()) return oldUnitDedx;

  if(G4CCIonChargeState::IonisationModelApplies(p, mat) &&
     !G4CCUseLowEnergyVanillaFallback(p, mat, kineticEnergy) &&
     G4CCCurrentChargeConfig::Instance().UseAtimaQDedxTable()) {
    const G4double eMeVPerU = G4CCIonChargeState::KineticEnergyPerNucleon(p, kineticEnergy)/CLHEP::MeV;
    const G4double w = G4CCCurrentChargeConfig::GlobalTransitionWeight(eMeVPerU);

    if(w > 0.999999) {
      G4double tabUnitMassDedx = 0.0;
      if(G4CCAtimaQDedxTable::Instance().GetMaterialUnitMassDEDX(p, mat, kineticEnergy,
                                                                  tabUnitMassDedx)) {
        return tabUnitMassDedx*(CLHEP::MeV*CLHEP::cm2/CLHEP::g)*mat->GetDensity();
      }
    }
  }

  return oldUnitDedx;
}


G4bool G4CCLindhardSorensenIonModel::ComputeDiracRestrictedDEDX(
  const G4Material* mat,const G4ParticleDefinition* p,G4double energy,
  G4double cut,G4double charge,G4double& dedx,G4double* fullOut,G4double* hardOut)
{
  dedx=0.;
  if(fullOut) *fullOut=0.;
  if(hardOut) *hardOut=0.;
  if(!G4CCCurrentChargeConfig::Instance().IsDiracMode() ||
     !G4CCIonChargeState::IonisationModelApplies(p,mat)) return false;
  const double eu=G4CCIonChargeState::KineticEnergyPerNucleon(p,energy)/MeV;
  if(eu<2.) return false;
  G4CCDiracPhysics::CheckCut(cut);
  SetParticle(p);
  const int z=G4CCIonChargeState::GetZp(p);
  const int n=G4CCIonChargeState::GetQe(p,charge);
  const double w=G4CCCurrentChargeConfig::GlobalTransitionWeight(eu);
  double tabFull=0.,tabHard=0.,fallbackFull=0.,fallbackHard=0.;
  if(w>0.) {
    tabFull=G4CCDiracPhysics::FullStopping(mat,p,n,energy);
    tabHard=G4CCDiracPhysics::HardStopping(mat,p,n,energy,cut);
  }
  if(w<1.) {
    fCorrections->GetParticleCharge(p,mat,energy);
    const double qeff2=fCorrections->EffectiveChargeSquareRatio(p,mat,energy);
    const bool high=eu>=G4CCCurrentChargeConfig::High0MeVPerU();
    if(!high) {
      if(eu<30.) {
        if(!G4CCAtimaEquilibriumDedxTable::Instance().GetMaterialDEDXPerVolume(p,mat,energy,fallbackFull)) {
          G4Exception("ComputeDiracRestrictedDEDX","CCDiracLowATIMA",FatalException,
                      "Required low-energy equilibrium ATIMA data are absent.");
          return false;
        }
      } else fallbackFull=G4CCDiracPhysics::EquilibriumStopping(mat,p,energy);
    } else {
      const double tau=energy/fMass, gamma=1.+tau;
      const double beta2=tau*(tau+2.)/(gamma*gamma);
      const double tmax=MaxSecondaryEnergy(p,energy), total=energy+fMass;
      const double prefactor=twopi_mc2_rcl2*qeff2*mat->GetElectronDensity()/beta2;
      const double base=ComputeOldUnitDEDXPerVolume(mat,p,energy,cut)*qeff2;
      const double bornCut=std::min(std::min(cut,tmax),fTlimit);
      double born=0.;
      if(bornCut>0. && bornCut<tmax) {
        born=std::log(tmax/bornCut)-beta2*(1.-bornCut/tmax);
        if(fSpin>0.) born+=(tmax*tmax-bornCut*bornCut)/(4.*total*total);
      }
      const double zDefinition=std::max(1.,std::abs(fDefinitionChargeInEplus));
      const double barkas=G4CCBarkasEnabled() ?
        2.*fCorrections->BarkasCorrection(p,mat,energy)*(zDefinition-1.)/zDefinition : 0.;
      const double zb=G4CCChooseCorrectionCharge(G4CCBlochZMode(),charge,z,fCorrections,p,mat,energy);
      const double bloch=G4CCDirectBlochEnabled() ? 2.*G4CCDirectBlochDeltaL(zb,std::sqrt(beta2)) : 0.;
      const int zm=G4CCChooseMottZForDeltaSubtraction(charge,z,fCorrections,p,mat,energy);
      const double mott=G4CCMottCorrection::FullStoppingCorrection(zm,std::sqrt(beta2));
      const double finiteSize=G4CCFiniteSizeLossIntegral(p,energy,zm);
      fallbackFull=base+prefactor*(born+barkas+bloch+mott-finiteSize);
    }
    fallbackHard=G4CCDeltaMeanLossPerVolumeForChargeState(mat,p,energy,cut,z,n,charge,
                 qeff2,0.,true,fCorrections,high);
  }
  const double full=w*tabFull+(1.-w)*fallbackFull;
  const double hard=w*tabHard+(1.-w)*fallbackHard;
  dedx=full-hard;
  if(!std::isfinite(dedx) || dedx < -1.e-12*std::max(full,hard)) {
    G4ExceptionDescription msg;
    msg<<"Negative/nonfinite restricted loss: Z="<<z<<", N="<<n<<", Eu="<<eu
       <<", full="<<full<<", hard="<<hard;
    G4Exception("ComputeDiracRestrictedDEDX","CCDiracBalance",FatalException,msg);
    return false;
  }
  dedx=std::max(0.,dedx);
  if(fullOut) *fullOut=full;
  if(hardOut) *hardOut=hard;
  return true;
}

G4bool G4CCLindhardSorensenIonModel::ComputeAtimaRestrictedDEDX(
  const G4Material* mat, const G4ParticleDefinition* p,
  G4double kineticEnergy, G4double deltaCut, G4double currentChargeInEplus,
  G4double& dedx, G4double* fullDedx, G4double* hardDedx, G4bool allowPlateau)
{
  if(G4CCCurrentChargeConfig::Instance().IsDiracMode())
    return ComputeDiracRestrictedDEDX(mat,p,kineticEnergy,deltaCut,currentChargeInEplus,
                                     dedx,fullDedx,hardDedx);
  dedx = 0.0;
  if(fullDedx != nullptr) { *fullDedx = 0.0; }
  if(hardDedx != nullptr) { *hardDedx = 0.0; }
  if(!G4CCCurrentChargeConfig::Instance().UseAtimaEquilibriumDedx() ||
     !G4CCIonChargeState::IonisationModelApplies(p, mat)) { return false; }
  const G4double eu = G4CCIonChargeState::KineticEnergyPerNucleon(p, kineticEnergy)/MeV;
  if(eu < 2.0 || (eu > 40.0 &&
     (!allowPlateau || eu > 1500.0 ||
      G4CCCurrentChargeConfig::GlobalTransitionWeight(eu) < 1.0))) { return false; }
  if(!(deltaCut > 0.0) || !std::isfinite(deltaCut) ||
     !std::isfinite(currentChargeInEplus)) {
    G4Exception("G4CCLindhardSorensenIonModel::ComputeAtimaRestrictedDEDX",
                "CCAtimaInvalidInput", FatalException,
                "ATIMA stopping requires a positive finite electron cut and finite charge.");
    return false;
  }
  SetParticle(p);
  const G4int zp = G4CCIonChargeState::GetZp(p);
  const G4int qe = std::max(0, std::min(zp,
    static_cast<G4int>(std::lround(zp - currentChargeInEplus))));
  fCorrections->GetParticleCharge(p, mat, kineticEnergy);
  const G4double qg2 = std::max(0.0, fCorrections->EffectiveChargeSquareRatio(p, mat, kineticEnergy));
  const G4double w = G4CCCurrentChargeConfig::GlobalTransitionWeight(eu);
  G4double eqFull = 0.0, eqHard = 0.0, ccFull = 0.0, ccHard = 0.0;
  if(w < 1.0) {
    if(!G4CCAtimaEquilibriumDedxTable::Instance().GetMaterialDEDXPerVolume(
         p, mat, kineticEnergy, eqFull)) {
      G4Exception("G4CCLindhardSorensenIonModel::ComputeAtimaRestrictedDEDX",
                  "CCAtimaMissingEquilibrium", FatalException,
                  "Missing equilibrium ATIMA stopping inside the configured 2--40 MeV/u domain.");
      return false;
    }
    eqHard = G4CCDeltaMeanLossPerVolumeForChargeState(mat, p, kineticEnergy,
      deltaCut, zp, qe, currentChargeInEplus, qg2, 0.0, true, fCorrections);
  }
  if(w > 0.0) {
    if(!G4CCCurrentChargeConfig::Instance().UseAtimaQDedxTable() ||
       !G4CCAtimaQDedxTable::Instance().GetMaterialDEDXPerVolume(p, mat,
         kineticEnergy, currentChargeInEplus*currentChargeInEplus, ccFull)) {
      G4Exception("G4CCLindhardSorensenIonModel::ComputeAtimaRestrictedDEDX",
                  "CCAtimaMissingCurrentCharge", FatalException,
                  "The 30--40 MeV/u transition requires the current-charge ATIMA table.");
      return false;
    }
    ccHard = G4CCDeltaMeanLossPerVolumeForChargeState(mat, p, kineticEnergy,
      deltaCut, zp, qe, currentChargeInEplus, qg2, 1.0, true, fCorrections);
  }
  const G4double full = (1.0-w)*eqFull + w*ccFull;
  const G4double hard = (1.0-w)*eqHard + w*ccHard;
  dedx = (1.0-w)*(eqFull-eqHard) + w*(ccFull-ccHard);
  if(!std::isfinite(dedx) || dedx < -1.e-12*std::max(full, hard)) {
    G4ExceptionDescription message;
    message << "Negative/nonfinite ATIMA restricted stopping: Z=" << zp
            << " E/A=" << eu << " MeV/u, q=" << currentChargeInEplus
            << ", material=" << mat->GetName() << ", cut=" << deltaCut/keV
            << " keV, full=" << full/(MeV/mm) << ", hard=" << hard/(MeV/mm)
            << " MeV/mm. The explicit spectrum exceeds the full stopping.";
    G4Exception("G4CCLindhardSorensenIonModel::ComputeAtimaRestrictedDEDX",
                "CCAtimaNegativeResidual", FatalException, message);
    return false;
  }
  dedx = std::max(0.0, dedx);
  if(fullDedx != nullptr) { *fullDedx = full; }
  if(hardDedx != nullptr) { *hardDedx = hard; }
  return true;
}

G4bool G4CCLindhardSorensenIonModel::ApplyAtimaAlongStep(
  const G4Material* mat, const G4ParticleDefinition* p,
  G4double preKinEnergy, G4double cut, G4double length, G4double& eloss,
  G4bool allowPlateau)
{
  SetParticle(p);
  const G4double preEu = G4CCIonChargeState::KineticEnergyPerNucleon(p, preKinEnergy)/MeV;
  if(G4CCCurrentChargeConfig::Instance().UseAtimaEquilibriumDedx() &&
     G4CCIonChargeState::IonisationModelApplies(p, mat) &&
     preEu >= 2.0 && (G4CCCurrentChargeConfig::Instance().IsDiracMode() || preEu <= 40.0 || (allowPlateau && preEu <= 1500.0)) && length > 0.0) {
    G4double deltaCut = cut;
    if(CurrentCouple() != nullptr && InTransportContext(p, mat)) {
      const auto* cuts = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(1);
      const G4int idx = CurrentCouple()->GetIndex();
      if(cuts != nullptr && idx >= 0 && static_cast<std::size_t>(idx) < cuts->size()) {
        deltaCut = (*cuts)[static_cast<std::size_t>(idx)];
      }
    }
    const G4double charge = MayUseCachedCharge(p, mat) ? fCachedChargeInEplus
      : static_cast<G4double>(G4CCIonChargeState::GetZp(p));
    G4double boundary = G4CCCurrentChargeConfig::Instance().IsDiracMode()
      ? 2.0/931.494*p->GetPDGMass() : 2.0*MeV*G4CCIonChargeState::GetA(p);
    while(G4CCIonChargeState::KineticEnergyPerNucleon(p, boundary)/MeV < 2.0) {
      boundary = std::nextafter(boundary, preKinEnergy);
    }
    G4double stopping = 0.0;
    const G4bool available = ComputeAtimaRestrictedDEDX(mat, p, preKinEnergy,
      deltaCut, charge, stopping, nullptr, nullptr, true);
    if(available) {
      const auto updateStopping = [&](G4double energy) {
        if(ComputeAtimaRestrictedDEDX(mat, p, energy, deltaCut, charge, stopping,
                                     nullptr, nullptr, true)) { return true; }
        G4Exception("G4CCLindhardSorensenIonModel::ApplyAtimaAlongStep",
          "CCAtimaStepEvaluation", FatalException,
          "Stopping is unavailable inside the ATIMA part of a transport step.");
        return false;
      };
      G4double loss = stopping*length;
      for(G4int iteration = 0; iteration < 3; ++iteration) {
        const G4double midpoint = std::max(boundary, preKinEnergy-0.5*loss);
        if(!updateStopping(midpoint)) { return false; }
        loss = stopping*length;
      }
      const G4double toBoundary = preKinEnergy-boundary;
      if(loss > toBoundary) {
        if(!updateStopping(0.5*(preKinEnergy+boundary))) { return false; }
        const G4double boundaryLength = stopping > 0.0 ? toBoundary/stopping : length;
        if(boundaryLength < length) {
          G4double remainingLength = length-std::max(0.0, boundaryLength);
          G4double legacyEnergy = std::nextafter(boundary, 0.0);
          while(G4CCIonChargeState::KineticEnergyPerNucleon(p, legacyEnergy)/MeV >= 2.0) {
            legacyEnergy = std::nextafter(legacyEnergy, 0.0);
          }
          const G4double nativeBoundary = fLowEnergyLimit/fProtonMassRatio;
          const auto legacySeed = [&](G4double energy) {
            SetParticle(p);
            const G4double scaledEnergy = energy*fProtonMassRatio;
            fCorrections->GetParticleCharge(G4GenericIon::GenericIon(), mat, scaledEnergy);
            const G4double reference = ComputeOldUnitDEDXPerVolume(
              mat, G4GenericIon::GenericIon(), scaledEnergy, cut);
            SetParticle(p);
            fCorrections->GetParticleCharge(p, mat, energy);
            return reference*GetChargeSquareRatio(p, mat, energy);
          };
          const auto legacyStepLoss = [&](G4double energy, G4double distance) {
            G4double result = legacySeed(energy)*distance;
            CorrectionsAlongStep(mat, p, energy, cut, distance, result);
            return result;
          };
          G4int substeps = 0;
          while(remainingLength > 0.0 && legacyEnergy > 0.0) {
            const G4double seed = legacySeed(legacyEnergy);
            if(!std::isfinite(seed) || seed <= 0.0 || ++substeps > 100000) {
              G4Exception("G4CCLindhardSorensenIonModel::ApplyAtimaAlongStep",
                "CCAtimaLegacyIntegration", FatalException,
                "Invalid or nonconvergent legacy integration below the ATIMA endpoint.");
              return false;
            }
            const G4double maximumLoss = 0.001*legacyEnergy;
            G4double distance = std::min(remainingLength, maximumLoss/seed);
            G4double stepLoss = legacyStepLoss(legacyEnergy, distance);
            for(G4int shrink = 0; stepLoss > maximumLoss && shrink < 8; ++shrink) {
              distance *= 0.99*maximumLoss/stepLoss;
              stepLoss = legacyStepLoss(legacyEnergy, distance);
            }
            if(!std::isfinite(stepLoss) || stepLoss < 0.0 || distance <= 0.0
               || stepLoss > maximumLoss*(1.0 + 1.e-12)) {
              G4Exception("G4CCLindhardSorensenIonModel::ApplyAtimaAlongStep",
                "CCAtimaLegacyStep", FatalException,
                "Invalid legacy substep below the ATIMA endpoint.");
              return false;
            }
            const G4double toNativeBoundary = legacyEnergy-nativeBoundary;
            if(toNativeBoundary > 0.0 && stepLoss >= toNativeBoundary) {
              G4double lowerDistance = 0.0, upperDistance = distance;
              for(G4int iteration = 0; iteration < 48; ++iteration) {
                const G4double middleDistance = 0.5*(lowerDistance+upperDistance);
                if(legacyStepLoss(legacyEnergy, middleDistance) < toNativeBoundary) {
                  lowerDistance = middleDistance;
                } else {
                  upperDistance = middleDistance;
                }
              }
              distance = 0.5*(lowerDistance+upperDistance);
              legacyEnergy = std::nextafter(nativeBoundary, 0.0);
            } else {
              legacyEnergy = std::max(0.0, legacyEnergy-stepLoss);
            }
            remainingLength = std::max(0.0, remainingLength-distance);
          }
          loss = preKinEnergy-legacyEnergy;
        }
      }
      eloss = std::max(0.0, std::min(preKinEnergy, loss));
      return true;
    }
  }
  return false;
}

void G4CCLindhardSorensenIonModel::CorrectionsAlongStep(
  const G4Material* mat, const G4ParticleDefinition* p,
  const G4double preKinEnergy, const G4double cut,
  const G4double& length, G4double& eloss)
{
  SetParticle(p);

  if(!InTransportContext(p, mat)) {
    G4double stopping = 0.0;
    if(ComputeAtimaRestrictedDEDX(mat, p, preKinEnergy, cut,
         static_cast<G4double>(G4CCIonChargeState::GetZp(p)), stopping,
         nullptr, nullptr, true)) {
      eloss = std::max(0.0, length)*stopping;
      return;
    }
  }

  const G4double preEu = G4CCIonChargeState::KineticEnergyPerNucleon(p, preKinEnergy)/MeV;
  if(ApplyAtimaAlongStep(mat, p, preKinEnergy, cut, length, eloss, false)) { return; }

  if(eloss >= preKinEnergy) { return; }

  const G4bool currentChargeApplies = G4CCIonChargeState::IonisationModelApplies(p, mat);
  const G4double eDensity = mat->GetElectronDensity();
  const G4double e = std::max(preKinEnergy - 0.5*eloss, 0.5*preKinEnergy);
  if(preEu > 40.0 &&
     G4CCIonChargeState::KineticEnergyPerNucleon(p, e)/MeV < 40.0 &&
     ApplyAtimaAlongStep(mat, p, preKinEnergy, cut, length, eloss, true)) { return; }
  const G4double tmax = MaxSecondaryEnergy(p, e);
  const G4double escaled = e*fProtonMassRatio;
  const G4double tau = e/fMass;

  const G4int zProjectile = G4CCIonChargeState::GetZp(p);
  const G4double eMeVPerU = G4CCIonChargeState::KineticEnergyPerNucleon(p, e)/CLHEP::MeV;
  const G4bool useLowEnergyVanillaFallback = currentChargeApplies && eMeVPerU <= 30.0;
  const G4bool inCurrentChargeDomain = currentChargeApplies && !useLowEnergyVanillaFallback;
  const G4double globalWeight = useLowEnergyVanillaFallback
    ? 0.0
    : G4CCCurrentChargeConfig::GlobalTransitionWeight(eMeVPerU);
  const G4bool mottActiveForStep = G4CCLSUseMottDeltaRuntime()
    && G4CCCurrentChargeConfig::MottEnabledForEnergy(eMeVPerU);

  G4double q2 = 0.0;
  G4double q2Current = 0.0;
  G4double q2Effective = 0.0;
  G4double qCurrentInEplus = 0.0;
  fCorrections->GetParticleCharge(p, mat, e);
  if(inCurrentChargeDomain) {
    qCurrentInEplus = MayUseCachedCharge(p, mat)
      ? fCachedChargeInEplus
      : static_cast<G4double>((zProjectile > 0) ? zProjectile : G4CCIonChargeState::GetZp(p));
    q2Current = std::max(0.0, qCurrentInEplus*qCurrentInEplus);
    q2Effective = std::max(0.0, fCorrections->EffectiveChargeSquareRatio(p, mat, e));
    q2 = std::max(0.0, globalWeight*q2Current + (1.0 - globalWeight)*q2Effective);
  } else {
    q2 = fCorrections->EffectiveChargeSquareRatio(p, mat, e);
    q2Current = q2;
    q2Effective = q2;
    qCurrentInEplus = std::sqrt(std::max(0.0, q2));
  }

  G4bool tableBlendActive = false;
  G4double tableEloss = 0.0;
  G4double tableFullDedx = 0.0;
  G4double tableHardDeltaDedx = 0.0;
  const auto checkResidual = [&](G4double residual, G4double full, G4double hard) {
    if(!std::isfinite(residual) || residual < -1.e-12*std::max(std::abs(full), std::abs(hard))) {
      G4ExceptionDescription message;
      message << "Negative/nonfinite Legacy restricted stopping: Z=" << zProjectile
              << ", E/A=" << eMeVPerU << " MeV/u, q=" << qCurrentInEplus
              << ", material=" << mat->GetName() << ", cut=" << cut/keV
              << " keV. The explicit spectrum exceeds the full stopping.";
      G4Exception("G4CCLindhardSorensenIonModel::CorrectionsAlongStep",
                  "CCAtimaNegativeResidual", FatalException, message);
      return false;
    }
    return true;
  };
  const G4int zMottLowTransition = inCurrentChargeDomain
    ? G4CCChooseMottZForDeltaSubtraction(qCurrentInEplus, zProjectile,
                                       fCorrections, p, mat, e) : 0;
  const G4bool correctLowTransition = inCurrentChargeDomain
    && eMeVPerU > 30.0 && eMeVPerU < 40.0 && mottActiveForStep
    && G4CCMottCorrection::IsAvailable(zMottLowTransition);
  const G4bool correctHighFallback = inCurrentChargeDomain
    && eMeVPerU >= 40.0 && globalWeight < 1.0 && mottActiveForStep
    && G4CCMottCorrection::IsAvailable(zMottLowTransition);
  const G4bool correctAnalyticFallback = correctLowTransition || correctHighFallback;

  G4double deltaCut = cut;
  if((correctAnalyticFallback ||
      (G4CCCurrentChargeConfig::Instance().UseAtimaEquilibriumDedx() &&
       inCurrentChargeDomain && G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(eMeVPerU))) &&
     CurrentCouple() != nullptr && InTransportContext(p, mat)) {
    const auto* cuts = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(1);
    const G4int idx = CurrentCouple()->GetIndex();
    if(cuts != nullptr && idx >= 0 && static_cast<std::size_t>(idx) < cuts->size()) {
      deltaCut = (*cuts)[static_cast<std::size_t>(idx)];
    }
  }
  if(inCurrentChargeDomain && globalWeight > 0.0 &&
     G4CCCurrentChargeConfig::Instance().UseAtimaQDedxTable()) {
    if(G4CCAtimaQDedxTable::Instance().GetMaterialDEDXPerVolume(p, mat, e, q2Current, tableFullDedx)) {
      const G4int qeCurrent = std::max(0, std::min(zProjectile,
        static_cast<G4int>(std::lround(static_cast<G4double>(zProjectile) - qCurrentInEplus))));

      tableHardDeltaDedx = G4CCDeltaMeanLossPerVolumeForChargeState(
        mat, p, e, deltaCut, zProjectile, qeCurrent, qCurrentInEplus, q2Effective,
        1.0, mottActiveForStep, fCorrections);
      const G4double tableResidual = tableFullDedx - tableHardDeltaDedx;
      const G4double tableContDedx = tableResidual;

      tableBlendActive = true;
      tableEloss = tableContDedx*length;

      if(preEu > 40.0 && preKinEnergy-tableEloss < 40.0*MeV*G4CCIonChargeState::GetA(p) &&
         ApplyAtimaAlongStep(mat, p, preKinEnergy, cut, length, eloss, true)) { return; }

      if(G4CCAtimaDeltaSubDiagnosticsEnabled()) {
        const G4int n = gAtimaDeltaSubDiagnosticPrints.fetch_add(1);
        if(n < 50) {
          G4cout << "G4CC_ATIMA_DELTA_SUB"
                 << " Zp=" << zProjectile
                 << " Qe=" << qeCurrent
                 << " q=" << qCurrentInEplus
                 << " EMeVPerU=" << eMeVPerU
                 << " cutMeV=" << cut/CLHEP::MeV
                 << " w=" << globalWeight
                 << " fullDedx(MeV/mm)=" << tableFullDedx/(CLHEP::MeV/CLHEP::mm)
                 << " hardDeltaDedx(MeV/mm)=" << tableHardDeltaDedx/(CLHEP::MeV/CLHEP::mm)
                 << " contDedx(MeV/mm)=" << tableContDedx/(CLHEP::MeV/CLHEP::mm)
                 << G4endl;
        }
      }

      if(globalWeight > 0.999999 && !correctAnalyticFallback) {
        if(!checkResidual(tableResidual, tableFullDedx, tableHardDeltaDedx)) { return; }
        eloss = std::max(0.0, std::min(tableEloss, preKinEnergy));
        return;
      }

      const G4double oldUnitDedx = ComputeOldUnitDEDXPerVolume(mat, p, e, cut);
      eloss = oldUnitDedx*q2Effective*length;
      q2 = q2Effective;
    }
  }

  if(correctAnalyticFallback) {
    const G4double gamma = tau + 1.0;
    const G4double beta2 = tau*(tau + 2.0)/(gamma*gamma);
    const G4double totalEnergy = e + fMass;
    const G4double prefactor = CLHEP::twopi_mc2_rcl2*q2*eDensity/beta2;
    const G4double restrictedBase = ComputeOldUnitDEDXPerVolume(mat, p, e, cut)*q2;
    const G4double bornCut = std::min(std::min(cut, tmax), fTlimit);
    G4double bornTailBracket = 0.0;
    if(bornCut > 0.0 && bornCut < tmax) {
      bornTailBracket = G4Log(tmax/bornCut) - beta2*(1.0 - bornCut/tmax);
      if(fSpin > 0.0) {
        bornTailBracket += (tmax*tmax - bornCut*bornCut)/(4.0*totalEnergy*totalEnergy);
      }
    }
    const G4double hardBornDedx = prefactor*bornTailBracket;
    const G4double chargeForBarkas = std::max(1.0, std::abs(fDefinitionChargeInEplus));
    const G4double barkas = G4CCBarkasEnabled()
      ? 2.0*fCorrections->BarkasCorrection(p, mat, e)*(chargeForBarkas - 1.0)/chargeForBarkas
      : 0.0;
    const G4double zBloch = G4CCChooseCorrectionCharge(G4CCBlochZMode(), qCurrentInEplus,
                                                      zProjectile, fCorrections, p, mat, e);
    const G4double twiceBloch = G4CCDirectBlochEnabled()
      ? 2.0*G4CCDirectBlochDeltaL(zBloch, std::sqrt(beta2)) : 0.0;
    const G4double fullMott = G4CCMottCorrection::FullStoppingCorrection(zMottLowTransition,
                                                                      std::sqrt(beta2));
    const G4int qe = std::max(0, std::min(zProjectile,
      static_cast<G4int>(std::lround(static_cast<G4double>(zProjectile) - qCurrentInEplus))));
    const G4double hardMottDedx = G4CCDeltaMeanLossPerVolumeForChargeState(
      mat, p, e, deltaCut, zProjectile, qe, qCurrentInEplus, q2Effective,
      tableBlendActive ? 0.0 : globalWeight, true, fCorrections,
      correctHighFallback);
    const G4double finiteSizeDedx = correctHighFallback
      ? CLHEP::twopi_mc2_rcl2*eDensity*q2Effective/beta2
        *(tableBlendActive ? 1.0 : 1.0-globalWeight)
        *G4CCFiniteSizeLossIntegral(p, e, zMottLowTransition) : 0.0;
    const G4double fullDedx = restrictedBase + hardBornDedx
      + prefactor*(barkas + twiceBloch + fullMott) - finiteSizeDedx;
    const G4double continuousDedx = fullDedx - hardMottDedx;
    G4double correctedLoss = continuousDedx*length;
    if(tableBlendActive) {
      correctedLoss = (1.0 - globalWeight)*correctedLoss + globalWeight*tableEloss;
    }
    const G4double fullLoss = length*(tableBlendActive
      ? (1.0-globalWeight)*fullDedx + globalWeight*tableFullDedx : fullDedx);
    const G4double hardLoss = length*(tableBlendActive
      ? (1.0-globalWeight)*hardMottDedx + globalWeight*tableHardDeltaDedx : hardMottDedx);
    if(!checkResidual(correctedLoss, fullLoss, hardLoss)) { return; }
    eloss = std::max(0.0, std::min(correctedLoss, preKinEnergy));
    if(G4CCBlochDiagnosticsEnabled()) {
      const G4int n = gBlochDiagnosticPrints.fetch_add(1);
      if(n < 50) {
        G4cout << (correctHighFallback ? "G4CC_UPPER_FALLBACK_STOPPING" : "G4CC_LOW_TRANSITION_STOPPING")
               << " EMeVPerU=" << eMeVPerU << " w=" << globalWeight
               << " atima=" << tableBlendActive << " q2=" << q2
               << " zBloch=" << zBloch << " zMott=" << zMottLowTransition
               << " cutBornMeV=" << bornCut/CLHEP::MeV
               << " cutDeltaMeV=" << deltaCut/CLHEP::MeV
               << " Barkas=" << barkas << " twiceBloch=" << twiceBloch
               << " fullMott=" << fullMott
               << " restrictedBase=" << restrictedBase/(CLHEP::MeV/CLHEP::mm)
               << " HB=" << hardBornDedx/(CLHEP::MeV/CLHEP::mm)
               << " HM=" << hardMottDedx/(CLHEP::MeV/CLHEP::mm)
               << " fullDedx=" << fullDedx/(CLHEP::MeV/CLHEP::mm)
               << " contDedx=" << continuousDedx/(CLHEP::MeV/CLHEP::mm)
               << G4endl;
      }
    }
    return;
  }

  const G4int zForData = (zProjectile > 0) ? zProjectile : p->GetAtomicNumber();

  G4double res = 0.0;
  if(escaled <= fLowEnergyLimit) {
    if(zForData > 2 && zForData <= 80) {
      res = fIonData->GetDEDX(mat, zForData, escaled, G4Log(escaled));
    }
    if(res > 0.0) {
      if(cut < tmax) {
        const G4double x = cut/tmax;
        res += (G4Log(x)*(tau + 1.0)*(tau + 1.0)/(tau*(tau + 2.0)) + 1.0 - x)
          *q2*CLHEP::twopi_mc2_rcl2*eDensity;
      }
      res *= length;
    } else {
      res = (fChargeSquareRatio > 0.0) ? eloss*q2/fChargeSquareRatio : eloss;
    }
  } else {
    const G4double gamma = tau + 1.0;
    const G4double beta2 = tau*(tau + 2.0)/(gamma*gamma);

    const G4int zForLS = (zProjectile > 0) ? zProjectile : fZin;

    const G4String zModeBloch = G4CCBlochZMode();
    const G4double zCorrBloch = G4CCChooseCorrectionCharge(zModeBloch, qCurrentInEplus,
                                                           zForLS, fCorrections, p, mat, e);

    const G4String zModeBarkas = G4CCBarkasZMode();
    const G4double zCorrBarkas = G4CCChooseCorrectionCharge(zModeBarkas, qCurrentInEplus,
                                                            zForLS, fCorrections, p, mat, e);

    const G4int zMottForContinuous = G4CCChooseMottZForDeltaSubtraction(qCurrentInEplus,
                                                                        zForLS,
                                                                        fCorrections,
                                                                        p, mat, e);

    G4double deltaL0 = 0.0;
    if(useLowEnergyVanillaFallback || !inCurrentChargeDomain || G4CCBarkasEnabled()) {
      const G4double chargeForBarkas = std::max(1.0, std::abs(fDefinitionChargeInEplus));
      deltaL0 = 2.0*fCorrections->BarkasCorrection(p, mat, e)
        *(chargeForBarkas - 1.0)/chargeForBarkas;
    }

    G4double deltaL = fLSData->GetDeltaL(zForLS, gamma);

    const G4bool suppressFullLS = inCurrentChargeDomain
      && mottActiveForStep
      && G4CCLSDisableContinuousLSRuntime()
      && G4CCMottCorrection::IsAvailable(zMottForContinuous);

    if(suppressFullLS) {
      if(G4CCLSDiagnosticsEnabled()) {
        const G4int n = gLSSuppressionDiagnosticPrints.fetch_add(1);
        if(n < 20) {
          G4cout << "G4CC_MOTT_DIAGNOSTICS continuous LS suppressed"
                 << " Zp=" << zForLS
                 << " zMottContinuous=" << zMottForContinuous
                 << " gamma=" << gamma
                 << " beta2=" << beta2
                 << " oldDeltaL=" << deltaL
                 << " keptDeltaL0=" << deltaL0
                 << G4endl;
        }
      }
      deltaL = 0.0;
    }

    G4double deltaLBloch = 0.0;
    if(inCurrentChargeDomain && mottActiveForStep && G4CCDirectBlochEnabled()) {
      deltaLBloch = G4CCDirectBlochDeltaL(zCorrBloch, std::sqrt(beta2));
    }

    G4double extraEnergy = 0.0;
    if(inCurrentChargeDomain && G4CCGeant4HighOrderFEnabled()) {
      extraEnergy += fCorrections->HighOrderCorrections(p, mat, e, cut)*length;
    }

    if(G4CCBlochDiagnosticsEnabled()) {
      const G4int n = gBlochDiagnosticPrints.fetch_add(1);
      if(n < 50) {
        G4cout << "G4CC_BLOCH_DIAGNOSTICS"
               << " Zp=" << zForLS
               << " zModeBloch=" << zModeBloch
               << " zCorrBloch=" << zCorrBloch
               << " zModeBarkas=" << zModeBarkas
               << " zCorrBarkas=" << zCorrBarkas
               << " barkasChargeVanilla=" << std::max(1.0, std::abs(fDefinitionChargeInEplus))
               << " zMottContinuous=" << zMottForContinuous
               << " q2=" << q2
               << " globalWeight=" << globalWeight
               << " eMeVPerU=" << eMeVPerU
               << " beta=" << std::sqrt(beta2)
               << " gamma=" << gamma
               << " deltaL_G4LS=" << deltaL
               << " deltaL_Barkas=" << deltaL0
               << " deltaL_BlochDirect=" << deltaLBloch
               << " dE_Bloch(MeV)="
               << (CLHEP::twopi_mc2_rcl2*q2*eDensity*deltaLBloch*length/beta2)/CLHEP::MeV
               << " extraHighOrderF(MeV)=" << extraEnergy/CLHEP::MeV
               << G4endl;
      }
    }

    res = eloss
      + CLHEP::twopi_mc2_rcl2*q2*eDensity*(deltaL + deltaL0 + deltaLBloch)*length/beta2
      + extraEnergy;
  }

  if(tableBlendActive) {
    const G4double fullLoss = globalWeight*tableFullDedx*length + (1.0-globalWeight)*res;
    const G4double hardLoss = globalWeight*tableHardDeltaDedx*length;
    res = globalWeight*tableEloss + (1.0 - globalWeight)*res;
    if(!checkResidual(res, fullLoss, hardLoss)) { return; }
    eloss = std::max(0.0, std::min(res, preKinEnergy));
    return;
  }

  if(res > preKinEnergy || 2.0*res < eloss) { res = eloss; }
  eloss = res;
}

void G4CCLindhardSorensenIonModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* secondaries,
  const G4MaterialCutsCouple* couple,
  const G4DynamicParticle* dp,
  G4double cut,
  G4double maxEnergy)
{
  const G4Material* material = (couple == nullptr) ? nullptr : couple->GetMaterial();
  if(dp != nullptr &&
     G4CCIonChargeState::IonisationModelApplies(dp->GetDefinition(), material) &&
     G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(
       G4CCIonChargeState::KineticEnergyPerNucleon(dp->GetDefinition(),
                                                 dp->GetKineticEnergy())/CLHEP::MeV)) {
    return;
  }

  G4double kineticEnergy = dp->GetKineticEnergy();
  const G4double tmax = MaxSecondaryEnergy(dp->GetDefinition(), kineticEnergy);
  const G4double minKinEnergy = std::min(cut, tmax);
  const G4double maxKinEnergy = std::min(maxEnergy, tmax);
  if(minKinEnergy >= maxKinEnergy) { return; }

  const G4double totalEnergy = kineticEnergy + fMass;
  const G4double totalEnergy2 = totalEnergy*totalEnergy;
  const G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*fMass)/totalEnergy2;

  G4double deltaKinEnergy = 0.0;
  G4double f = 0.0;
  G4double f1 = 0.0;
  G4double fmax = 1.0;
  if(fSpin > 0.0) { fmax += 0.5*maxKinEnergy*maxKinEnergy/totalEnergy2; }

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  G4double rndm[2];

  do {
    rndmEngine->flatArray(2, rndm);
    deltaKinEnergy = minKinEnergy*maxKinEnergy
      /(minKinEnergy*(1.0 - rndm[0]) + maxKinEnergy*rndm[0]);

    f = 1.0 - beta2*deltaKinEnergy/tmax;
    if(fSpin > 0.0) {
      f1 = 0.5*deltaKinEnergy*deltaKinEnergy/totalEnergy2;
      f += f1;
    }
  } while(fmax*rndm[1] > f);

  const G4double x = fFormFactor*deltaKinEnergy;
  if(x > 1.e-6) {
    const G4double x1 = 1.0 + x;
    G4double grej = 1.0/(x1*x1);
    if(fSpin > 0.0) {
      const G4double x2 = 0.5*electron_mass_c2*deltaKinEnergy/(fMass*fMass);
      grej *= (1.0 + fMagMoment2*(x2 - f1/f)/(1.0 + x2));
    }
    if(rndmEngine->flat() > grej) { return; }
  }

  G4ThreeVector deltaDirection;
  if(UseAngularGeneratorFlag()) {
    const G4Material* mat = couple->GetMaterial();
    const G4int zTarget = SelectRandomAtomNumber(mat);
    deltaDirection = GetAngularDistribution()->SampleDirection(dp, deltaKinEnergy, zTarget, mat);
  } else {
    const G4double deltaMomentum =
      std::sqrt(deltaKinEnergy*(deltaKinEnergy + 2.0*electron_mass_c2));
    G4double cost = deltaKinEnergy*(totalEnergy + electron_mass_c2)
      /(deltaMomentum*dp->GetTotalMomentum());
    cost = std::min(cost, 1.0);
    const G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
    const G4double phi = CLHEP::twopi*rndmEngine->flat();
    deltaDirection.set(sint*std::cos(phi), sint*std::sin(phi), cost);
    deltaDirection.rotateUz(dp->GetMomentumDirection());
  }

  auto* delta = new G4DynamicParticle(fElectron, deltaDirection, deltaKinEnergy);
  secondaries->push_back(delta);

  kineticEnergy -= deltaKinEnergy;
  G4ThreeVector finalP = dp->GetMomentum() - delta->GetMomentum();
  finalP = finalP.unit();

  fParticleChangeForLoss->SetProposedKineticEnergy(kineticEnergy);
  fParticleChangeForLoss->SetProposedMomentumDirection(finalP);
}

G4double G4CCLindhardSorensenIonModel::MaxSecondaryEnergy(
  const G4ParticleDefinition* particleDefinition, G4double kineticEnergy)
{
  SetParticle(particleDefinition);
  const G4double tau = kineticEnergy/fMass;
  return 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.0)
    /(1.0 + 2.0*(tau + 1.0)*fElectronMassRatio
      + fElectronMassRatio*fElectronMassRatio);
}
