#ifndef G4CCDiracPhysics_hh
#define G4CCDiracPhysics_hh 1
#include "globals.hh"
#include <array>
class G4CCDiracTable;
class G4ParticleDefinition;
class G4Material;
namespace G4CCDiracPhysics {
  G4CCDiracTable& Tables();
  G4bool MaterialInDomain(const G4Material*);
  void RequireMaterialCoverage(G4int projectileZ,const G4Material*);
  G4bool Applies(const G4ParticleDefinition*, const G4Material*);
  void CheckCut(G4double cut);
  std::array<G4double,3> Moments(const G4ParticleDefinition*, G4int n,
                               G4double kineticEnergy, G4double cut);
  G4double Sample(const G4ParticleDefinition*, G4int n, G4double kineticEnergy,
                  G4double cut, G4double uNode, G4double uCDF);
  G4double FullStopping(const G4Material*, const G4ParticleDefinition*,
                        G4int n, G4double kineticEnergy);
  G4double EquilibriumStopping(const G4Material*, const G4ParticleDefinition*,
                               G4double kineticEnergy);
  G4double HardStopping(const G4Material*, const G4ParticleDefinition*, G4int n,
                        G4double kineticEnergy, G4double cut);
}
#endif
