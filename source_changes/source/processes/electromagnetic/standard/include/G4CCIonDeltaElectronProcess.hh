#ifndef G4CCIonDeltaElectronProcess_hh
#define G4CCIonDeltaElectronProcess_hh 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleChangeForLoss.hh"
#include <vector>

class G4VEmAngularDistribution;

class G4CCIonDeltaElectronProcess : public G4VDiscreteProcess
{
 public:
  explicit G4CCIonDeltaElectronProcess(const G4String& name = "CCIonDeltaElectron");
  ~G4CCIonDeltaElectronProcess() override;

  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  G4bool IsApplicable(const G4ParticleDefinition& particle) override;

  G4double GetMeanFreePath(const G4Track& track,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  G4VParticleChange* PostStepDoIt(const G4Track& track,
                                  const G4Step& step) override;

  static G4double FormFactorAcceptance(const G4ParticleDefinition* particle,
                                       G4double kineticEnergy,
                                       G4double deltaEnergy);
  static G4double FormFactorThresholdEnergy(const G4ParticleDefinition* particle,
                                            G4double kineticEnergy);

 private:
  static G4double MaxSecondaryEnergy(const G4ParticleDefinition* particle,
                                     G4double kineticEnergy);
  static G4double Shape(G4double deltaEnergy,
                        G4double beta2,
                        G4double tMax,
                        G4double totalEnergy2,
                        G4double spin);
  static G4double NumericalIntegral(G4int Zp,
                                    G4int zMott,
                                    G4int Qe,
                                    G4double tMin,
                                    G4double tMax,
                                    G4double beta2,
                                    G4double totalEnergy2,
                                    G4double spin,
                                    G4double qEffective2,
                                    G4double globalWeight,
                                    G4bool useMott);
  static G4double ElectronProductionCut(const G4MaterialCutsCouple* couple);

  G4double CrossSectionPerVolume(const G4Track& track) const;
  G4double CrossSectionPerVolume(const G4Track& track, G4double kineticEnergy,
                                 G4bool analyticOnly) const;
  G4double SampleDeltaEnergy(const G4Track& track,
                             G4double tMin,
                             G4double tMax,
                             G4double beta2,
                             G4double totalEnergy2,
                             G4double spin,
                             G4double qEffective2,
                             G4double globalWeight,
                             G4bool useMott) const;

  G4double AnalyticCrossSectionUpperBound(const G4Track&, G4double low, G4double high);
  G4double CandidateCrossSection(const G4Track&);
  struct Envelope {
    const G4ParticleDefinition* particle;
    const G4MaterialCutsCouple* couple;
    G4double cut;
    G4int charge, electrons;
    G4double maximum;
    std::vector<G4double> energies, bounds;
  };
  std::vector<Envelope> fEnvelopes;
  G4double fCandidateCrossSection=0.;
  G4ParticleChangeForLoss fParticleChange;
  G4VEmAngularDistribution* fAngularDistribution = nullptr;
};

#endif
