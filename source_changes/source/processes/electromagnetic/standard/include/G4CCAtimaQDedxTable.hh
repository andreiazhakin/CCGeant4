#ifndef G4CCAtimaQDedxTable_hh
#define G4CCAtimaQDedxTable_hh 1

#include "globals.hh"

#include <atomic>
#include <map>
#include <vector>

class G4Material;
class G4ParticleDefinition;

class G4CCAtimaQDedxTable
{
 public:
  static G4CCAtimaQDedxTable& Instance();

  static void CheckConfigurationChangeAllowed();

  const G4String& Path() const;
  void SetPath(const G4String& path);

  G4bool IsLoaded() const;
  G4bool IsAvailable() const;

  G4bool GetUnitMassDEDX(G4int zProjectile,
                         G4int zTarget,
                         G4double kineticEnergyMeVPerU,
                         G4double& valueMeVcm2PerG) const;

  G4bool GetMaterialUnitMassDEDX(const G4ParticleDefinition* particle,
                                 const G4Material* material,
                                 G4double kineticEnergy,
                                 G4double& valueMeVcm2PerG) const;

  G4bool GetMaterialDEDXPerVolume(const G4ParticleDefinition* particle,
                                  const G4Material* material,
                                  G4double kineticEnergy,
                                  G4double chargeSquare,
                                  G4double& dedxPerVolume) const;

 private:
  G4CCAtimaQDedxTable();

  struct Point
  {
    G4double eMeVPerU = 0.0;
    G4double unitMassDedx = 0.0;
  };

  using Curve = std::vector<Point>;

  void Clear();
  G4bool EnsureLoaded() const;
  G4bool Load() const;
  static long Key(G4int zProjectile, G4int zTarget);

  G4String fPath;

  mutable std::atomic<G4bool> fLoadAttempted{false};
  mutable std::atomic<G4bool> fLoaded{false};
  mutable std::map<long, Curve> fData;
};

#endif
