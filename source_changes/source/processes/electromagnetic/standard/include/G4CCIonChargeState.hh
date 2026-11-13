#ifndef G4CCIonChargeState_hh
#define G4CCIonChargeState_hh 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

class G4DynamicParticle;
class G4Material;
class G4ParticleDefinition;
class G4ParticleChangeForLoss;
class G4Track;

class G4CCIonChargeState
{
 public:
  static G4bool IsIonDefinition(const G4ParticleDefinition* particle);
  static G4bool IsIonTrack(const G4Track& track);

  static G4int GetZp(const G4ParticleDefinition* particle);
  static G4int GetA(const G4ParticleDefinition* particle);
  static G4int GetZp(const G4Track& track);
  static G4int GetA(const G4Track& track);

  static G4double GetChargeInEplus(const G4DynamicParticle* dynamicParticle);
  static G4double GetChargeInEplus(const G4Track& track);

  static G4int GetQe(const G4Track& track);
  static G4int GetQe(const G4ParticleDefinition* particle, G4double chargeInEplus);
  static G4int ClampQe(G4int Zp, G4int Qe);

  static G4double ChargeFromQe(G4int Zp, G4int Qe);
  static void ProposeQe(G4ParticleChangeForLoss& particleChange, G4int Zp, G4int Qe);
  static void ProposeBareIon(G4ParticleChangeForLoss& particleChange, G4int Zp);

  static G4double KineticEnergyPerNucleon(const G4Track& track);
  static G4double KineticEnergyPerNucleon(const G4ParticleDefinition* particle,
                                          G4double kineticEnergy);

  static G4int MaxTargetZ();
  static G4bool MaterialIsSupportedTarget(const G4Material* material, G4int maxZt = 96);
  static G4bool MaterialIsLightOnly(const G4Material* material, G4int maxZt = 96);
  static G4bool CurrentChargeModelApplies(const G4Track& track);
  static G4bool CurrentChargeModelApplies(const G4ParticleDefinition* particle,
                                          const G4Material* material);

  static G4bool IonisationModelApplies(const G4Track& track);
  static G4bool IonisationModelApplies(const G4ParticleDefinition*,const G4Material*);

  static G4double CurrentChargeSquare(const G4Track& track);

  static G4double EffectiveChargeSquare(const G4Track& track);

  static G4double StoppingChargeSquare(const G4Track& track);

  static G4double ChargeSquareRatioToBareZ(const G4Track& track);
};

#endif
