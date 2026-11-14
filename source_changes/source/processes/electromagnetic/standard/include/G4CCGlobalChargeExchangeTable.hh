#ifndef G4CCGlobalChargeExchangeTable_hh
#define G4CCGlobalChargeExchangeTable_hh 1

#include "globals.hh"

#include <atomic>
#include <cstdint>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class G4CCGlobalChargeExchangeTable
{
 public:
  struct Transition
  {
    G4int Zt = 0;
    G4int QeIn = 0;
    G4int QeOut = 0;
    G4int dQe = 0;
    G4double sigmaCm2 = 0.0;
  };

  static G4CCGlobalChargeExchangeTable& Instance();

  void SetPath(const G4String& path);
  const G4String& Path() const;

  void EnsureLoaded();
  G4bool IsLoaded() const;

  void RequirePair(G4int Zp, G4int Zt);
  G4bool HasPair(G4int Zp, G4int Zt) const;

  std::vector<Transition> FindTransitions(G4int Zp,
                                           G4int A,
                                           G4int Zt,
                                           G4double kineticEnergyMeVPerU,
                                           G4int QeIn);

  void FindTransitions(G4int Zp,
                       G4int A,
                       G4int Zt,
                       G4double kineticEnergyMeVPerU,
                       G4int QeIn,
                       std::vector<Transition>& out);

  G4double InterpolatedSigmaCm2(G4int Zp,
                                G4int A,
                                G4int Zt,
                                G4int QeIn,
                                G4int dQe,
                                G4double kineticEnergyMeVPerU);

  std::size_t NumberOfCurves() const;
  std::size_t NumberOfPoints() const;

 private:
  G4CCGlobalChargeExchangeTable() = default;

  struct Point
  {
    G4float energyMeVPerU = 0.0F;
    G4int QeOut = 0;
    G4float sigmaCm2 = 0.0F;
  };

  struct Curve
  {
    G4int dQe = 0;
    G4int QeOut = 0;
    std::vector<Point> points;
  };

  static std::uint32_t PackGroupKey(G4int Zp, G4int Zt, G4int QeIn);

  static G4double Interpolate(const std::vector<Point>& points,
                              G4double kineticEnergyMeVPerU);
  static G4String Trim(const G4String& text);

  void LoadUnlocked();
  void ValidatePairUnlocked(G4int Zp, G4int Zt) const;
  void CheckExplicitQuery(G4int Zp, G4int Zt, G4int QeIn,
                          G4double energyMeVPerU) const;

  mutable std::mutex fLoadMutex;
  G4String fPath;
  std::atomic_bool fLoaded{false};
  std::atomic_bool fLoadAttempted{false};
  G4bool fStrict = false;
  std::unordered_set<std::uint32_t> fCompletePairs;
  std::unordered_map<std::uint32_t, std::vector<Curve>> fData;
  std::atomic_size_t fNumberOfCurves{0};
  std::atomic_size_t fNumberOfPoints{0};
};

#endif
