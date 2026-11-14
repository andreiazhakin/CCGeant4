#ifndef G4CCDiracTable_hh
#define G4CCDiracTable_hh 1

#include "globals.hh"
#include <array>
#include <atomic>
#include <memory>
#include <mutex>

class G4CCDiracTable
{
 public:
  static G4CCDiracTable& Instance();
  ~G4CCDiracTable();
  void Configure(const G4String& deltaPath, const G4String& stoppingPath);
  void Validate() const;
  G4bool HasProjectile(G4int z) const;
  G4bool HasPair(G4int z, G4int zt) const;
  G4double EnergyMin() const;
  G4double EnergyMax() const;
  static constexpr G4double MinTransfer() { return 0.00099; }
  static constexpr G4double TabulatedMinTransfer() { return 0.001; }
  static G4double MaxTransfer(G4double energyMeVPerU);

  std::array<G4double, 3> DeltaMoments(G4int z, G4int n, G4double energyMeVPerU,
                                      G4double cutMeV, G4double maximumTransferMeV = 0.) const;
  G4double DeltaCrossSectionUpperBound(G4int z, G4int n,
    G4double energyLow, G4double energyHigh, G4double cutMeV,
    G4double maximumTransferLowMeV, G4double maximumTransferHighMeV) const;
  G4double DeltaDensity(G4int z, G4int n, G4double energyMeVPerU,
                        G4double transferMeV, G4double maximumTransferMeV = 0.) const;
  G4double SampleDelta(G4int z, G4int n, G4double energyMeVPerU,
                       G4double cutMeV, G4double uNode, G4double uCDF,
                       G4double maximumTransferMeV = 0.) const;
  G4double FullStopping(G4int z, G4int zt, G4int n, G4double energyMeVPerU) const;
  G4double EquilibriumStopping(G4int z, G4int zt, G4double energyMeVPerU,
                               G4bool leftLimit = false) const;
  G4double SoftCoefficient(G4int z, G4int zt, G4double energyMeVPerU,
                            G4bool leftLimit = false) const;
  std::size_t MemoryBytes() const;

 private:
  G4CCDiracTable();
  struct Impl;
  mutable std::unique_ptr<Impl> fData;
  mutable std::mutex fMutex;
  mutable std::atomic<G4bool> fLoaded{false};
  G4String fDeltaPath, fStoppingPath;
};

#endif
