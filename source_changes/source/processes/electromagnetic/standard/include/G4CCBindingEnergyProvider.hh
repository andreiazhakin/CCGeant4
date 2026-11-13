#ifndef G4CCBindingEnergyProvider_hh
#define G4CCBindingEnergyProvider_hh 1

#include "globals.hh"
#include <memory>

class G4CCBindingEnergyProvider
{
 public:
  virtual ~G4CCBindingEnergyProvider() = default;

  virtual G4double BindingEnergy(G4int Zp, G4int electronIndex) const = 0;
};

class G4CCHydrogenicBindingEnergyProvider final : public G4CCBindingEnergyProvider
{
 public:
  G4double BindingEnergy(G4int Zp, G4int electronIndex) const override;
};

class G4CCBindingEnergyProviderStore
{
 public:
  static G4CCBindingEnergyProviderStore& Instance();

  const G4CCBindingEnergyProvider& Provider() const;
  void SetProvider(std::unique_ptr<G4CCBindingEnergyProvider> provider);

 private:
  G4CCBindingEnergyProviderStore();
  std::unique_ptr<G4CCBindingEnergyProvider> fProvider;
};

#endif
