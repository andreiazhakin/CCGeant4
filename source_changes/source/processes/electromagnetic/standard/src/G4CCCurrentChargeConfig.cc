#include "G4CCCurrentChargeConfig.hh"
#include "G4CCDiracTable.hh"

#include "G4CCAtimaQDedxTable.hh"
#include "G4CCGlobalChargeExchangeTable.hh"
#include "G4CCScreening.hh"
#include "G4Exception.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <memory>
#include <string>
#include <filesystem>
#include <stdexcept>


namespace
{
std::unique_ptr<G4CCPhysicsOptions> gPendingOptions;
G4CCCurrentChargeConfig* gExistingConfig = nullptr;

[[noreturn]] void G4CCConfigFail(const G4String& message)
{
  G4Exception("G4CCCurrentChargeConfig", "CCExplicitConfig001", FatalException, message.c_str());
  throw std::runtime_error(message);
}

void G4CCCheckPreInit()
{
  if(G4Threading::IsWorkerThread() ||
     G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit) {
    G4CCConfigFail("Explicit physics options must be set on the master in PreInit.");
  }
}

void G4CCValidateOptions(const G4CCPhysicsOptions& o)
{
  const auto& t = o.transition;
  if(!std::isfinite(t.low0) || !std::isfinite(t.low1) ||
     !std::isfinite(t.high0) || !std::isfinite(t.high1) ||
     !(30*MeV <= t.low0 && t.low0 < t.low1 && t.low1 <= t.high0 &&
       t.high0 < t.high1 && t.high1 <= 2000*MeV)) {
    G4CCConfigFail("Require 30 <= LOW0 < LOW1 <= HIGH0 < HIGH1 <= 2000 MeV/u (speed coordinate).");
  }
  if(o.tables.delta.empty() || o.tables.stopping.empty() ||
     o.tables.atimaLow.empty() || o.tables.global.empty()) {
    G4CCConfigFail("Dirac physics requires explicit delta, stopping, atimaLow and GLOBAL CSV paths.");
  }
  if(o.scattering.mode != G4CCScatteringMode::Opt4MultipleScattering &&
     o.scattering.mode != G4CCScatteringMode::NativeSingleScattering &&
     o.scattering.mode != G4CCScatteringMode::TabulatedSingleScattering) {
    G4CCConfigFail("Invalid ion scattering mode.");
  }
  if(o.fluctuationMode != "mix" && o.fluctuationMode != "ion" &&
     o.fluctuationMode != "universal") {
    G4CCConfigFail("fluctuationMode must be mix, ion or universal.");
  }
  if(!std::isfinite(o.fluctuationMixPower) || !(o.fluctuationMixPower > 0.)) {
    G4CCConfigFail("fluctuationMixPower must be finite and positive.");
  }
}

G4String G4CCConfigLower(G4String value)
{
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

G4bool G4CCConfigEnvFlagValue(const char* value, G4bool defaultValue)
{
  if(value == nullptr) { return defaultValue; }
  const G4String s = G4CCConfigLower(G4String(value));
  if(s == "0" || s == "false" || s == "off" || s == "no") { return false; }
  if(s == "1" || s == "true" || s == "on" || s == "yes") { return true; }
  return defaultValue;
}

G4bool G4CCConfigEnvFlag(const char* name, G4bool defaultValue)
{
  return G4CCConfigEnvFlagValue(std::getenv(name), defaultValue);
}

G4double G4CCConfigEnvDouble(const char* name, G4double defaultValue)
{
  const char* value = std::getenv(name);
  if(value == nullptr || value[0] == '\0') { return defaultValue; }
  try {
    return std::stod(value);
  } catch(...) {
    return defaultValue;
  }
}

G4double G4CCConfigClamp01(G4double x)
{
  return std::max(0.0, std::min(1.0, x));
}

G4double G4CCConfigSmoothStep(G4double x)
{
  x = G4CCConfigClamp01(x);
  return x*x*(3.0 - 2.0*x);
}

G4double G4CCConfigLogRamp(G4double eMeVPerU, G4double e0, G4double e1)
{
  if(eMeVPerU <= 0.0) { return 0.0; }
  if(e1 <= e0 || e0 <= 0.0) { return (eMeVPerU >= e0) ? 1.0 : 0.0; }
  const G4double x = (std::log(eMeVPerU) - std::log(e0))/(std::log(e1) - std::log(e0));
  return G4CCConfigSmoothStep(x);
}

struct G4CCGlobalWeightRuntimeConfig
{
  G4bool enableSmoothFallback = true;
  G4double low0 = 30.0;
  G4double low1 = 40.0;
  G4double high0 = 1500.0;
  G4double high1 = 2000.0;
};

const G4CCGlobalWeightRuntimeConfig& G4CCGlobalWeightRuntime()
{
  static const G4CCGlobalWeightRuntimeConfig cfg = []() {
    G4CCGlobalWeightRuntimeConfig c;
    c.enableSmoothFallback = G4CCConfigEnvFlag("G4CC_ENABLE_SMOOTH_GLOBAL_FALLBACK", true);
    c.low0 = G4CCConfigEnvDouble("G4CC_GLOBAL_WEIGHT_LOW0_MEV_U", 30.0);
    c.low1 = G4CCConfigEnvDouble("G4CC_GLOBAL_WEIGHT_LOW1_MEV_U", 40.0);
    c.high0 = G4CCConfigEnvDouble("G4CC_GLOBAL_WEIGHT_HIGH0_MEV_U", 1500.0);
    c.high1 = G4CCConfigEnvDouble("G4CC_GLOBAL_WEIGHT_HIGH1_MEV_U", 2000.0);
    return c;
  }();
  return cfg;
}
}

G4CCCurrentChargeConfig& G4CCCurrentChargeConfig::Instance()
{
  static G4CCCurrentChargeConfig config;
  return config;
}

G4CCCurrentChargeConfig::G4CCCurrentChargeConfig()
{
  if(gPendingOptions) {
    SetOptions(*gPendingOptions);
    gPendingOptions.reset();
    gExistingConfig = this;
    return;
  }
  const char* env = std::getenv("G4CC_GLOBAL_TABLE");
  if(env != nullptr) { fTablePath = env; }

  const char* atimaDedx = std::getenv("G4CC_ATIMA_Q_DEDX_TABLE");
  if(atimaDedx == nullptr) { atimaDedx = std::getenv("G4CC_ATIMA_DEDX_TABLE"); }
  if(atimaDedx != nullptr) { fAtimaQDedxTablePath = atimaDedx; }

  fUseMottDeltaElectrons = G4CCConfigEnvFlag("G4CC_USE_MOTT_DELTA", true);

  const G4double mottMinimum = G4CCConfigEnvDouble("G4CC_MOTT_MIN_ENERGY_MEV_U",
                                                 fUseAtimaEquilibriumDedx ? 2.0 : 30.0);
  if(fUseAtimaEquilibriumDedx &&
     (!fUseMottDeltaElectrons || !std::isfinite(mottMinimum) || mottMinimum > 2.0)) {
    G4Exception("G4CCCurrentChargeConfig::G4CCCurrentChargeConfig", "CCConfig001",
                FatalException,
                "Equilibrium ATIMA stopping requires G4CC_USE_MOTT_DELTA=1 and "
                "G4CC_MOTT_MIN_ENERGY_MEV_U <= 2.");
  }

  const auto& globalWeight = G4CCGlobalWeightRuntime();
  if(fUseAtimaEquilibriumDedx &&
     (!globalWeight.enableSmoothFallback || globalWeight.low0 != 30.0 || globalWeight.low1 != 40.0)) {
    G4Exception("G4CCCurrentChargeConfig::G4CCCurrentChargeConfig", "CCConfig003",
                FatalException,
                "Equilibrium ATIMA stopping requires smooth GLOBAL weighting and "
                "G4CC_GLOBAL_WEIGHT_LOW0_MEV_U=30, G4CC_GLOBAL_WEIGHT_LOW1_MEV_U=40.");
  }

  const char* disableLS = std::getenv("G4CC_DISABLE_CONTINUOUS_LS_FOR_MOTT_DELTA");
  if(disableLS != nullptr) {
    const G4String value(disableLS);
    fDisableContinuousLSWhenMottDelta = !(value == "0" || value == "false" || value == "FALSE" || value == "off" || value == "OFF");
  }

  const char* fluctMode = std::getenv("G4CC_FLUCTUATION_MODE");
  if(fluctMode != nullptr && fluctMode[0] != '\0') {
    SetFluctuationMode(G4String(fluctMode));
  }
  const char* mixPower = std::getenv("G4CC_FLUCTUATION_MIX_POWER");
  if(mixPower != nullptr && mixPower[0] != '\0') {
    SetFluctuationMixPower(G4CCConfigEnvDouble("G4CC_FLUCTUATION_MIX_POWER", 1.0));
  }
  gExistingConfig = this;
}

void G4CCCurrentChargeConfig::Configure(const G4CCPhysicsOptions& options)
{
  G4CCCheckPreInit();
  G4CCValidateOptions(options);
  if(gExistingConfig) {
    gExistingConfig->SetOptions(options);
  } else {
    gPendingOptions = std::make_unique<G4CCPhysicsOptions>(options);
    (void)Instance();
  }
}

void G4CCCurrentChargeConfig::SetOptions(const G4CCPhysicsOptions& options)
{
  G4CCCheckPreInit();
  if(fFrozen.load(std::memory_order_acquire)) {
    G4CCConfigFail("Physics options are frozen after ConstructProcess; create a new run configuration.");
  }
  G4CCValidateOptions(options);
  fOptions = options;
  fDiracMode = true;
  G4CCDiracTable::Instance().Configure(options.tables.delta,options.tables.stopping);
  fTablePath = options.tables.global;
  G4CCGlobalChargeExchangeTable::Instance().SetPath(fTablePath);
  fAtimaQDedxTablePath.clear();
  fUseAtimaQDedxTable = false;
  fUseAtimaEquilibriumDedx = true;
  fUseMottDeltaElectrons = true;
  fDisableContinuousLSWhenMottDelta = true;
  fFluctuationMode = options.fluctuationMode;
  fFluctuationMixPower = options.fluctuationMixPower;
}

void G4CCCurrentChargeConfig::Freeze()
{
  if(!fDiracMode || fFrozen.load(std::memory_order_acquire)) { return; }
  G4CCValidateOptions(fOptions);
  for(const auto& path : {fOptions.tables.delta, fOptions.tables.stopping,
                          fOptions.tables.atimaLow, fOptions.tables.global}) {
    if(!std::filesystem::exists(std::filesystem::path(path.c_str()))) {
      G4CCConfigFail("Required table path does not exist: " + path);
    }
  }
  fFrozen.store(true, std::memory_order_release);
}

void G4CCCurrentChargeConfig::CheckLegacyMutation(const char* method) const
{
  if(fDiracMode) {
    G4CCConfigFail(G4String(method) + ": use constructor options for explicit Dirac physics; legacy setters cannot override them.");
  }
}

G4double G4CCCurrentChargeConfig::Low0MeVPerU()
{
  const auto& c = Instance();
  return c.IsDiracMode() ? c.Options().transition.low0/MeV : G4CCGlobalWeightRuntime().low0;
}
G4double G4CCCurrentChargeConfig::Low1MeVPerU()
{
  const auto& c = Instance();
  return c.IsDiracMode() ? c.Options().transition.low1/MeV : G4CCGlobalWeightRuntime().low1;
}
G4double G4CCCurrentChargeConfig::High0MeVPerU()
{
  const auto& c = Instance();
  return c.IsDiracMode() ? c.Options().transition.high0/MeV : G4CCGlobalWeightRuntime().high0;
}
G4double G4CCCurrentChargeConfig::High1MeVPerU()
{
  const auto& c = Instance();
  return c.IsDiracMode() ? c.Options().transition.high1/MeV : G4CCGlobalWeightRuntime().high1;
}

const G4String& G4CCCurrentChargeConfig::TablePath() const
{
  return fTablePath;
}

void G4CCCurrentChargeConfig::SetTablePath(const G4String& path)
{
  CheckLegacyMutation("SetTablePath");
  fTablePath = path;
  G4CCGlobalChargeExchangeTable::Instance().SetPath(path);
}

const G4String& G4CCCurrentChargeConfig::AtimaQDedxTablePath() const
{
  return fAtimaQDedxTablePath;
}

void G4CCCurrentChargeConfig::SetAtimaQDedxTablePath(const G4String& path)
{
  CheckLegacyMutation("SetAtimaQDedxTablePath");
  G4CCAtimaQDedxTable::Instance().SetPath(path);
  fAtimaQDedxTablePath = path;
}

G4bool G4CCCurrentChargeConfig::UseAtimaQDedxTable() const
{
  return fUseAtimaQDedxTable;
}

G4bool G4CCCurrentChargeConfig::UseAtimaEquilibriumDedx() const
{
  return fUseAtimaEquilibriumDedx;
}

G4bool G4CCCurrentChargeConfig::UseCCDeltaAtEnergy(G4double kineticEnergyMeVPerU)
{
  if(!std::isfinite(kineticEnergyMeVPerU)) { return false; }
  if(Instance().UseAtimaEquilibriumDedx()) { return kineticEnergyMeVPerU >= 2.0; }
  static const G4double legacyMaximum =
    std::max(0.0, G4CCConfigEnvDouble("G4CC_LOW_ENERGY_VANILLA_MAX_MEV_U", 30.0));
  return kineticEnergyMeVPerU > legacyMaximum;
}

G4double G4CCCurrentChargeConfig::ScreeningExponent() const
{
  return G4CCScreening::Exponent();
}

void G4CCCurrentChargeConfig::SetScreeningExponent(G4double exponent)
{
  CheckLegacyMutation("SetScreeningExponent");
  G4CCScreening::SetExponent(exponent);
}

G4bool G4CCCurrentChargeConfig::UseMottDeltaElectrons() const
{
  return fUseMottDeltaElectrons;
}

void G4CCCurrentChargeConfig::SetUseMottDeltaElectrons(G4bool value)
{
  CheckLegacyMutation("SetUseMottDeltaElectrons");
  if(fUseAtimaEquilibriumDedx && !value) {
    G4Exception("G4CCCurrentChargeConfig::SetUseMottDeltaElectrons", "CCConfig002",
                FatalException,
                "Cannot disable Mott delta electrons in equilibrium ATIMA mode.");
    return;
  }
  fUseMottDeltaElectrons = value;
}

G4bool G4CCCurrentChargeConfig::DisableContinuousLSWhenMottDelta() const
{
  return fDisableContinuousLSWhenMottDelta;
}

void G4CCCurrentChargeConfig::SetDisableContinuousLSWhenMottDelta(G4bool value)
{
  CheckLegacyMutation("SetDisableContinuousLSWhenMottDelta");
  fDisableContinuousLSWhenMottDelta = value;
}

const G4String& G4CCCurrentChargeConfig::FluctuationMode() const
{
  return fFluctuationMode;
}

void G4CCCurrentChargeConfig::SetFluctuationMode(const G4String& mode)
{
  CheckLegacyMutation("SetFluctuationMode");
  G4String m = G4CCConfigLower(mode);
  if(m == "mixture" || m == "mixed" || m == "auto") { m = "mix"; }
  if(m != "universal" && m != "ion" && m != "mix") { m = "universal"; }
  fFluctuationMode = m;
}

G4double G4CCCurrentChargeConfig::FluctuationMixPower() const
{
  return fFluctuationMixPower;
}

void G4CCCurrentChargeConfig::SetFluctuationMixPower(G4double power)
{
  CheckLegacyMutation("SetFluctuationMixPower");
  if(!(power > 0.0) || !std::isfinite(power)) { power = 1.0; }
  fFluctuationMixPower = std::max(0.05, std::min(16.0, power));
}

G4double G4CCCurrentChargeConfig::IonFluctuationMixProbability(G4double kineticEnergyMeVPerU)
{
  const G4double w = GlobalTransitionWeight(kineticEnergyMeVPerU);
  const G4double pBase = G4CCConfigClamp01(1.0 - w);
  const G4double power = Instance().FluctuationMixPower();
  return G4CCConfigClamp01(std::pow(pBase, power));
}

G4double G4CCCurrentChargeConfig::GlobalTransitionWeight(G4double kineticEnergyMeVPerU)
{
  if(Instance().IsDiracMode()) {
    if(!std::isfinite(kineticEnergyMeVPerU)) { return 0.; }
    const G4double lower = G4CCConfigLogRamp(kineticEnergyMeVPerU, Low0MeVPerU(), Low1MeVPerU());
    const G4double upper = 1. - G4CCConfigLogRamp(kineticEnergyMeVPerU, High0MeVPerU(), High1MeVPerU());
    return G4CCConfigClamp01(lower*upper);
  }
  if(!std::isfinite(kineticEnergyMeVPerU) || kineticEnergyMeVPerU <= 30.0) { return 0.0; }
  const auto& cfg = G4CCGlobalWeightRuntime();
  if(!cfg.enableSmoothFallback) { return 1.0; }

  const G4double lower = G4CCConfigLogRamp(kineticEnergyMeVPerU, cfg.low0, cfg.low1);
  const G4double upper = 1.0 - G4CCConfigLogRamp(kineticEnergyMeVPerU, cfg.high0, cfg.high1);
  return G4CCConfigClamp01(lower*upper);
}

G4double G4CCCurrentChargeConfig::MottMinimumEnergyMeVPerU()
{
  if(Instance().IsDiracMode()) { return 2.; }
  static const G4double threshold =
    std::max(0.0, G4CCConfigEnvDouble("G4CC_MOTT_MIN_ENERGY_MEV_U",
                                   Instance().UseAtimaEquilibriumDedx() ? 2.0 : 30.0));
  return threshold;
}

G4bool G4CCCurrentChargeConfig::MottEnabledForEnergy(G4double kineticEnergyMeVPerU)
{
  const G4double threshold = MottMinimumEnergyMeVPerU();
  return threshold <= 0.0 || kineticEnergyMeVPerU >= threshold;
}

void G4CCCurrentChargeConfig::EnsureMessenger()
{
  if(!fMessenger) { fMessenger = std::make_unique<G4CCCurrentChargeMessenger>(this); }
}

G4CCCurrentChargeMessenger::G4CCCurrentChargeMessenger(G4CCCurrentChargeConfig* config)
  : fConfig(config)
{
  fDirectory = new G4UIdirectory("/currentCharge/");
  fDirectory->SetGuidance("Current-charge heavy-ion model controls.");

  fTableCmd = new G4UIcmdWithAString("/currentCharge/table", this);
  fTableCmd->SetGuidance("Set CSV file with GLOBAL charge-exchange cross sections.");
  fTableCmd->SetParameterName("path", false);
  fTableCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fAtimaQDedxTableCmd = new G4UIcmdWithAString("/currentCharge/atimaQDedxTable", this);
  fAtimaQDedxTableCmd->SetGuidance("Set CSV file with ATIMA/GLOBAL q-scaled unit-q2 stopping table.");
  fAtimaQDedxTableCmd->SetParameterName("path", false);
  fAtimaQDedxTableCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fAtimaQDedxTableCmd->SetToBeBroadcasted(false);

  fExponentCmd = new G4UIcmdWithADouble("/currentCharge/screeningExponent", this);
  fExponentCmd->SetGuidance("Set exponent a in f(x)=1/(1+x^a) for delta-electron screening.");
  fExponentCmd->SetParameterName("a", false);
  fExponentCmd->SetDefaultValue(2.0);
  fExponentCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fUseMottDeltaCmd = new G4UIcmdWithABool("/currentCharge/useMottDeltaElectrons", this);
  fUseMottDeltaCmd->SetGuidance("Use Mott/Rutherford correction for explicit heavy-ion delta-electron production in the current-charge domain.");
  fUseMottDeltaCmd->SetParameterName("use", false);
  fUseMottDeltaCmd->SetDefaultValue(true);
  fUseMottDeltaCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDisableLSCmd = new G4UIcmdWithABool("/currentCharge/disableContinuousLSWhenMottDelta", this);
  fDisableLSCmd->SetGuidance("Suppress the full continuous Lindhard-Sorensen correction in the current-charge domain when Mott delta-electron production is active and Mott data are available.");
  fDisableLSCmd->SetParameterName("disable", false);
  fDisableLSCmd->SetDefaultValue(true);
  fDisableLSCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fFluctuationModeCmd = new G4UIcmdWithAString("/currentCharge/fluctuationMode", this);
  fFluctuationModeCmd->SetGuidance("Energy-loss fluctuation fallback mode: universal, ion, or mix.");
  fFluctuationModeCmd->SetParameterName("mode", false);
  fFluctuationModeCmd->SetDefaultValue("mix");
  fFluctuationModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fFluctuationMixPowerCmd = new G4UIcmdWithADouble("/currentCharge/fluctuationMixPower", this);
  fFluctuationMixPowerCmd->SetGuidance("Power in P(IonFluctuations)=(1-w)^power for fluctuationMode=mix.");
  fFluctuationMixPowerCmd->SetParameterName("power", false);
  fFluctuationMixPowerCmd->SetDefaultValue(1.0);
  fFluctuationMixPowerCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

G4CCCurrentChargeMessenger::~G4CCCurrentChargeMessenger()
{
  delete fFluctuationMixPowerCmd;
  delete fFluctuationModeCmd;
  delete fDisableLSCmd;
  delete fUseMottDeltaCmd;
  delete fExponentCmd;
  delete fAtimaQDedxTableCmd;
  delete fTableCmd;
  delete fDirectory;
}

void G4CCCurrentChargeMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
  if(command == fTableCmd) {
    fConfig->SetTablePath(value);
  } else if(command == fAtimaQDedxTableCmd) {
    fConfig->SetAtimaQDedxTablePath(value);
  } else if(command == fExponentCmd) {
    fConfig->SetScreeningExponent(fExponentCmd->GetNewDoubleValue(value));
  } else if(command == fUseMottDeltaCmd) {
    fConfig->SetUseMottDeltaElectrons(fUseMottDeltaCmd->GetNewBoolValue(value));
  } else if(command == fDisableLSCmd) {
    fConfig->SetDisableContinuousLSWhenMottDelta(fDisableLSCmd->GetNewBoolValue(value));
  } else if(command == fFluctuationModeCmd) {
    fConfig->SetFluctuationMode(value);
  } else if(command == fFluctuationMixPowerCmd) {
    fConfig->SetFluctuationMixPower(fFluctuationMixPowerCmd->GetNewDoubleValue(value));
  }
}
