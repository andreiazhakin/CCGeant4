#ifndef G4CC_TABULATED_ION_SCATTERING_HH
#define G4CC_TABULATED_ION_SCATTERING_HH

#include "G4VEmProcess.hh"
#include "G4CCElasticTableData.hh"
#include <memory>

class G4CCTabulatedIonScattering final : public G4VEmProcess {
 public:
  explicit G4CCTabulatedIonScattering(
      std::shared_ptr<const G4CCElasticTableData> data = {},
      G4bool useTables = true, G4double recoilThreshold = 0.,
      G4bool allowEnvironment = true);
  ~G4CCTabulatedIonScattering() override = default;

  G4bool IsApplicable(const G4ParticleDefinition&) override;
  void PreparePhysicsTable(const G4ParticleDefinition&) override;
  void StartTracking(G4Track*) override;
  G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                           const G4Material*) override;
  void ProcessDescription(std::ostream&) const override;

 protected:
  void InitialiseProcess(const G4ParticleDefinition*) override;
  void StreamProcessInfo(std::ostream&) const override;

 private:
  G4bool useTables_;
  G4bool isInitialised_ = false;
};

#endif
