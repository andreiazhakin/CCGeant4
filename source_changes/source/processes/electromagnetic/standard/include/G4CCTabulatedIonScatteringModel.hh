#ifndef G4CC_TABULATED_ION_SCATTERING_MODEL_HH
#define G4CC_TABULATED_ION_SCATTERING_MODEL_HH

#include "G4VEmModel.hh"
#include "G4CCElasticTableData.hh"
#include "G4ThreeVector.hh"
#include <memory>

class G4Track;

namespace tabulated {


struct ElasticKinematics {
  double finalKineticEnergy = 0;
  double recoilKineticEnergy = 0;
  double momentumCM = 0;
  double q2 = 0;
  G4ThreeVector primaryMomentum;
  G4ThreeVector recoilMomentum;
};

ElasticKinematics ComputeElasticKinematics(double projectileMass,
                                           double targetMass,
                                           double kineticEnergy,
                                           double x,
                                           double phi);

class G4CCTabulatedIonScatteringModel final : public G4VEmModel {
 public:
  explicit G4CCTabulatedIonScatteringModel(
      std::shared_ptr<const ElasticTableData> data, bool useTables = true);
  ~G4CCTabulatedIonScatteringModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;
  void InitialiseLocal(const G4ParticleDefinition*, G4VEmModel*) override;
  void StartTracking(G4Track*) override;

  G4double CrossSectionPerVolume(const G4Material*, const G4ParticleDefinition*,
                                G4double kineticEnergy,
                                G4double cutEnergy = 0,
                                G4double maxEnergy = DBL_MAX) override;
  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                     G4double kineticEnergy, G4double Z,
                                     G4double A = 0,
                                     G4double cutEnergy = 0,
                                     G4double maxEnergy = DBL_MAX) override;
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                         const G4MaterialCutsCouple*,
                         const G4DynamicParticle*, G4double = 0,
                         G4double = DBL_MAX) override;
  void ModelDescription(std::ostream&) const override;

  void SetRecoilThreshold(G4double);
  G4double CrossSectionForIsotope(const G4ParticleDefinition*, G4double energy,
                                 int Z, int A, int isomerLevel = 0);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

using TabulatedIonElasticModel = G4CCTabulatedIonScatteringModel;
}

using G4CCTabulatedIonScatteringModel =
    tabulated::G4CCTabulatedIonScatteringModel;
#endif
