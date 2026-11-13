#ifndef G4CCLindhardSorensenIonModel_h
#define G4CCLindhardSorensenIonModel_h 1

#include <vector>

#include "G4VEmModel.hh"
#include "G4NistManager.hh"

class G4EmCorrections;
class G4ParticleChangeForLoss;
class G4LindhardSorensenData;
class G4BraggModel;
class G4BetheBlochModel;
class G4IonICRU73Data;
class G4Track;
class G4CCIonIonisation;

class G4CCLindhardSorensenIonModel : public G4VEmModel
{
 public:
  explicit G4CCLindhardSorensenIonModel(const G4ParticleDefinition* p = nullptr,
                                        const G4String& nam = "CCLindhardSorensen");

  ~G4CCLindhardSorensenIonModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double ChargeSquareRatio(const G4Track&) override;

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple* couple) override;

  G4double ComputeCrossSectionPerElectron(const G4ParticleDefinition*,
                                          G4double kineticEnergy,
                                          G4double cutEnergy,
                                          G4double maxEnergy);

  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                      G4double kineticEnergy,
                                      G4double Z,
                                      G4double A,
                                      G4double cutEnergy,
                                      G4double maxEnergy) override;

  G4double CrossSectionPerVolume(const G4Material*,
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy,
                                 G4double cutEnergy,
                                 G4double maxEnergy) override;

  G4double ComputeDEDXPerVolume(const G4Material*,
                                const G4ParticleDefinition*,
                                G4double kineticEnergy,
                                G4double cutEnergy) override;


  G4double GetChargeSquareRatio(const G4ParticleDefinition*,
                                const G4Material*,
                                G4double kineticEnergy) override;

  G4double GetParticleCharge(const G4ParticleDefinition*,
                             const G4Material*,
                             G4double kineticEnergy) override;

  void CorrectionsAlongStep(const G4Material*,
                            const G4ParticleDefinition*,
                            const G4double kineticEnergy,
                            const G4double cutEnergy,
                            const G4double& length,
                            G4double& eloss) override;

  G4bool ComputeDiracRestrictedDEDX(const G4Material*,const G4ParticleDefinition*,
    G4double kineticEnergy,G4double deltaCut,G4double currentChargeInEplus,
    G4double& dedx,G4double* fullDedx=nullptr,G4double* hardDedx=nullptr);

  G4bool ComputeAtimaRestrictedDEDX(const G4Material*,
                                  const G4ParticleDefinition*,
                                  G4double kineticEnergy,
                                  G4double deltaCut,
                                  G4double currentChargeInEplus,
                                  G4double& dedx,
                                  G4double* fullDedx = nullptr,
                                  G4double* hardDedx = nullptr,
                                  G4bool allowPlateau = false);

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                         const G4MaterialCutsCouple*,
                         const G4DynamicParticle*,
                         G4double tmin,
                         G4double maxEnergy) override;

  G4CCLindhardSorensenIonModel& operator=(const G4CCLindhardSorensenIonModel&) = delete;
  G4CCLindhardSorensenIonModel(const G4CCLindhardSorensenIonModel&) = delete;

 protected:
  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                              G4double kineticEnergy) override;

 private:
  friend class G4CCIonIonisation;
  class TransportScope {
   public:
    explicit TransportScope(const G4Track&);
    ~TransportScope();
    TransportScope(const TransportScope&) = delete;
    TransportScope& operator=(const TransportScope&) = delete;
   private:
    const G4Track* fPrevious = nullptr;
  };
  static G4bool InTransportContext(const G4ParticleDefinition*, const G4Material*);
  G4bool MayUseCachedCharge(const G4ParticleDefinition*, const G4Material*) const;

  G4bool ApplyAtimaAlongStep(const G4Material*, const G4ParticleDefinition*,
                            G4double kineticEnergy, G4double cutEnergy,
                            G4double length, G4double& eloss,
                            G4bool allowPlateau);
  void SetupParameters();
  inline void SetParticle(const G4ParticleDefinition*);
  G4double ComputeOldUnitDEDXPerVolume(const G4Material*,
                                       const G4ParticleDefinition*,
                                       G4double kineticEnergy,
                                       G4double cutEnergy);

  static G4IonICRU73Data* fIonData;
  static G4LindhardSorensenData* fLSData;

  const G4ParticleDefinition* fParticle = nullptr;
  G4ParticleDefinition* fElectron = nullptr;
  G4EmCorrections* fCorrections = nullptr;
  G4ParticleChangeForLoss* fParticleChangeForLoss = nullptr;
  G4NistManager* fNist = nullptr;
  G4BraggModel* fBraggModel = nullptr;
  G4BetheBlochModel* fBetheBlochModel = nullptr;

  G4int fZin = 1;
  G4double fMass = 0.0;
  G4double fTlimit = DBL_MAX;
  G4double fSpin = 0.0;
  G4double fMagMoment2 = 0.0;
  G4double fChargeSquareRatio = 1.0;
  G4double fDefinitionChargeInEplus = 1.0;
  G4double fElectronMassRatio = 0.0;
  G4double fProtonMassRatio = 1.0;
  G4double fFormFactor = 0.0;
  G4double fTwoLog10 = 0.0;
  G4double fLowEnergyLimit = 0.0;

  G4double fCachedChargeInEplus = 0.0;
  G4bool fHasCachedCurrentCharge = false;

  G4bool fOwnsSharedData = false;
};

inline void G4CCLindhardSorensenIonModel::SetParticle(const G4ParticleDefinition* p)
{
  if(fParticle != p) {
    fParticle = p;
    SetupParameters();
  }
}

#endif
