#ifndef G4CCCurrentChargeConfig_hh
#define G4CCCurrentChargeConfig_hh 1

#include "globals.hh"
#include "G4CCPhysicsOptions.hh"
#include "G4UImessenger.hh"

#include <memory>
#include <atomic>

class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIdirectory;
class G4UIcommand;

class G4CCCurrentChargeMessenger;

class G4CCCurrentChargeConfig
{
 public:
  static G4CCCurrentChargeConfig& Instance();

  static void Configure(const G4CCPhysicsOptions& options);
  void SetOptions(const G4CCPhysicsOptions& options);
  G4bool IsDiracMode() const { return fDiracMode; }
  const G4CCPhysicsOptions& Options() const { return fOptions; }
  void Freeze();
  static G4double Low0MeVPerU();
  static G4double Low1MeVPerU();
  static G4double High0MeVPerU();
  static G4double High1MeVPerU();

  const G4String& TablePath() const;
  void SetTablePath(const G4String& path);

  const G4String& AtimaQDedxTablePath() const;
  void SetAtimaQDedxTablePath(const G4String& path);

  G4bool UseAtimaQDedxTable() const;

  G4bool UseAtimaEquilibriumDedx() const;

  static G4bool UseCCDeltaAtEnergy(G4double kineticEnergyMeVPerU);

  G4double ScreeningExponent() const;
  void SetScreeningExponent(G4double exponent);

  G4bool UseMottDeltaElectrons() const;
  void SetUseMottDeltaElectrons(G4bool value);

  G4bool DisableContinuousLSWhenMottDelta() const;
  void SetDisableContinuousLSWhenMottDelta(G4bool value);

  const G4String& FluctuationMode() const;
  void SetFluctuationMode(const G4String& mode);
  G4double FluctuationMixPower() const;
  void SetFluctuationMixPower(G4double power);
  static G4double IonFluctuationMixProbability(G4double kineticEnergyMeVPerU);

  static G4double GlobalTransitionWeight(G4double kineticEnergyMeVPerU);

  static G4double MottMinimumEnergyMeVPerU();
  static G4bool MottEnabledForEnergy(G4double kineticEnergyMeVPerU);

  void EnsureMessenger();

 private:
  G4CCCurrentChargeConfig();
  void CheckLegacyMutation(const char* method) const;
  G4CCPhysicsOptions fOptions;
  G4bool fDiracMode = false;
  std::atomic_bool fFrozen{false};
  G4String fTablePath;
  G4String fAtimaQDedxTablePath;
  G4bool fUseAtimaQDedxTable = true;
  G4bool fUseAtimaEquilibriumDedx = true;
  G4bool fUseMottDeltaElectrons = true;
  G4bool fDisableContinuousLSWhenMottDelta = true;
  G4String fFluctuationMode = "mix";
  G4double fFluctuationMixPower = 1.0;
  std::unique_ptr<G4CCCurrentChargeMessenger> fMessenger;
};

class G4CCCurrentChargeMessenger : public G4UImessenger
{
 public:
  explicit G4CCCurrentChargeMessenger(G4CCCurrentChargeConfig* config);
  ~G4CCCurrentChargeMessenger() override;

  void SetNewValue(G4UIcommand* command, G4String value) override;

 private:
  G4CCCurrentChargeConfig* fConfig = nullptr;
  G4UIdirectory* fDirectory = nullptr;
  G4UIcmdWithAString* fTableCmd = nullptr;
  G4UIcmdWithAString* fAtimaQDedxTableCmd = nullptr;
  G4UIcmdWithADouble* fExponentCmd = nullptr;
  G4UIcmdWithABool* fUseMottDeltaCmd = nullptr;
  G4UIcmdWithABool* fDisableLSCmd = nullptr;
  G4UIcmdWithAString* fFluctuationModeCmd = nullptr;
  G4UIcmdWithADouble* fFluctuationMixPowerCmd = nullptr;
};

#endif
