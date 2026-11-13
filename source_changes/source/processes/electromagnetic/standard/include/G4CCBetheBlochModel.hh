#ifndef G4CCBetheBlochModel_hh
#define G4CCBetheBlochModel_hh 1

#include "G4BetheBlochModel.hh"

class G4CCBetheBlochModel : public G4BetheBlochModel
{
 public:
  explicit G4CCBetheBlochModel(const G4ParticleDefinition* p = nullptr,
                               const G4String& name = "CCBetheBloch");
  ~G4CCBetheBlochModel() override = default;

  G4double ChargeSquareRatio(const G4Track& track) override;

  G4double CrossSectionPerVolume(const G4Material* material,
                                 const G4ParticleDefinition* particle,
                                 G4double kineticEnergy,
                                 G4double cutEnergy,
                                 G4double maxEnergy) override;

  G4double ComputeDEDXPerVolume(const G4Material* material,
                                const G4ParticleDefinition* particle,
                                G4double kineticEnergy,
                                G4double cutEnergy) override;

  G4double GetChargeSquareRatio(const G4ParticleDefinition* particle,
                                const G4Material* material,
                                G4double kineticEnergy) override;

  G4double GetParticleCharge(const G4ParticleDefinition* particle,
                             const G4Material* material,
                             G4double kineticEnergy) override;

  void CorrectionsAlongStep(const G4Material* material,
                            const G4ParticleDefinition* particle,
                            const G4double kineticEnergy,
                            const G4double cutEnergy,
                            const G4double& length,
                            G4double& eloss) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>* secondaries,
                         const G4MaterialCutsCouple* couple,
                         const G4DynamicParticle* dynamicParticle,
                         G4double tmin,
                         G4double maxEnergy) override;

 private:
  G4double fCachedChargeInEplus = 0.0;
  G4bool fHasCachedCurrentCharge = false;
};

#endif
