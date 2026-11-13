#ifndef G4CCMixedFluctuations_hh
#define G4CCMixedFluctuations_hh 1

#include "G4VEmFluctuationModel.hh"

#include <memory>

class G4IonFluctuations;
class G4Material;
class G4UniversalFluctuation;

class G4CCMixedFluctuations : public G4VEmFluctuationModel
{
 public:
  explicit G4CCMixedFluctuations(const G4String& name = "CCMixedFluc");
  ~G4CCMixedFluctuations() override;

  G4double SampleFluctuations(const G4MaterialCutsCouple*,
                              const G4DynamicParticle*,
                              const G4double tcut,
                              const G4double tmax,
                              const G4double length,
                              const G4double meanLoss) override;

  G4double Dispersion(const G4Material*,
                      const G4DynamicParticle*,
                      const G4double tcut,
                      const G4double tmax,
                      const G4double length) override;

  void InitialiseMe(const G4ParticleDefinition*) override;
  void SetParticleAndCharge(const G4ParticleDefinition*, G4double q2) override;

  G4CCMixedFluctuations& operator=(const G4CCMixedFluctuations&) = delete;
  G4CCMixedFluctuations(const G4CCMixedFluctuations&) = delete;

 private:
  G4double IonProbability(const G4DynamicParticle*, const G4Material*) const;

  std::unique_ptr<G4UniversalFluctuation> fUniversal;
  std::unique_ptr<G4IonFluctuations> fIon;
};

#endif
