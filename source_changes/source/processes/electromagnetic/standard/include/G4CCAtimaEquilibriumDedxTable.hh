#ifndef G4CCAtimaEquilibriumDedxTable_hh
#define G4CCAtimaEquilibriumDedxTable_hh 1

#include "globals.hh"

#include <cstddef>
#include <map>
#include <vector>

class G4Material;
class G4ParticleDefinition;

class G4CCAtimaEquilibriumDedxTable
{
 public:
  static const G4CCAtimaEquilibriumDedxTable& Instance();

  const G4String& Path() const { return fPath; }
  std::size_t NumberOfCurves() const { return fData.size(); }
  std::size_t NumberOfRows() const;
  G4bool HasPair(G4int zProjectile, G4int zTarget) const;
  void RequirePair(G4int zProjectile, G4int zTarget) const;
  G4double GetReferenceAtomicMass(G4int zProjectile, G4int zTarget) const;

  G4bool GetMassDEDX(G4int zProjectile, G4int zTarget,
                    G4double energyMeVPerU, G4double& massDedx) const;

  G4bool GetMassDEDXLeftLimit(G4int zProjectile, G4int zTarget,
                             G4double energyMeVPerU, G4double& massDedx) const;

  G4bool GetMaterialDEDXPerVolume(const G4ParticleDefinition* particle,
                                 const G4Material* material,
                                 G4double kineticEnergy,
                                 G4double& dedxPerVolume) const;

 private:
  struct Point
  {
    G4double energy;
    G4double main;
    G4double left;
  };
  struct Curve
  {
    G4double referenceAtomicMass = 0.0;
    std::vector<Point> points;
  };
  using Data = std::map<G4int, Curve>;

  G4CCAtimaEquilibriumDedxTable();
  G4CCAtimaEquilibriumDedxTable(const G4CCAtimaEquilibriumDedxTable&) = delete;
  G4CCAtimaEquilibriumDedxTable& operator=(const G4CCAtimaEquilibriumDedxTable&) = delete;

  static G4String SelectPath();
  static Data Load(const G4String& path);
  static G4int Key(G4int zp, G4int zt) { return 1000*zp + zt; }
  G4bool Lookup(G4int zp, G4int zt, G4double energy,
                G4bool leftLimit, G4double& value) const;

  const G4String fPath;
  const Data fData;
};

#endif
