#include "G4EmStandardPhysics_option4_CurrentCharge.hh"

#include "G4CCCurrentChargeConfig.hh"
#include "G4CCNuclearStopping.hh"
#include "G4NuclearStopping.hh"
#include "G4EmParameters.hh"
#include "G4Exception.hh"
#include <exception>
#include "G4CCTabulatedIonScattering.hh"
#include "G4VMultipleScattering.hh"

#include <set>
#include <map>
#include "G4CCIonChargeExchangeProcess.hh"
#include "G4CCIonChargeState.hh"
#include "G4CCIonDeltaElectronProcess.hh"
#include "G4CCIonIonisation.hh"
#include "G4CCNonEMChargeResetProcess.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsConstructorFactory.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmStandardPhysics_option4_CurrentCharge);

namespace
{
G4bool HasProcessNamed(const G4ProcessManager* pm, const G4String& name)
{
  if(pm == nullptr || pm->GetProcessList() == nullptr) { return false; }
  const auto* list = pm->GetProcessList();
  for(std::size_t i = 0; i < list->size(); ++i) {
    const auto* p = (*list)[static_cast<G4int>(i)];
    if(p != nullptr && p->GetProcessName() == name) { return true; }
  }
  return false;
}

void ReplaceNuclearStoppingWithEnergyBudget()
{
  std::set<G4ProcessManager*> seen;
  std::map<G4VProcess*, G4CCNuclearStopping*> replacements;
  auto* iterator = G4ParticleTable::GetParticleTable()->GetIterator();
  iterator->reset();
  while((*iterator)()) {
    auto* pm = iterator->value()->GetProcessManager();
    if(!pm || !seen.insert(pm).second) continue;
    for(G4int i = pm->GetProcessListLength()-1; i >= 0; --i) {
      auto* old = dynamic_cast<G4NuclearStopping*>((*pm->GetProcessList())[i]);
      if(!old) continue;
      const G4int atRest = pm->GetProcessOrdering(old, idxAtRest);
      const G4int along = pm->GetProcessOrdering(old, idxAlongStep);
      const G4int post = pm->GetProcessOrdering(old, idxPostStep);
      auto& replacement = replacements[old];
      if(!replacement) {
        replacement = new G4CCNuclearStopping(old->GetProcessName());
        replacement->SetMaxKinEnergy(G4EmParameters::Instance()->MaxNIELEnergy());
        replacement->SetVerboseLevel(old->GetVerboseLevel());
      }
      pm->RemoveProcess(old);
      pm->AddProcess(replacement, atRest, along, post);
    }
  }
}

void RemoveProcessesNamed(G4ProcessManager* pm, const G4String& name)
{
  if(pm == nullptr || pm->GetProcessList() == nullptr) { return; }
  auto* list = pm->GetProcessList();
  for(G4int i = static_cast<G4int>(list->size()) - 1; i >= 0; --i) {
    auto* process = (*list)[i];
    if(process != nullptr && process->GetProcessName() == name) {
      auto* removed = pm->RemoveProcess(i);
      delete removed;
    }
  }
}
}

G4EmStandardPhysics_option4_CurrentCharge::G4EmStandardPhysics_option4_CurrentCharge(
  const G4CCPhysicsOptions& options, G4int ver, const G4String& name)
  : G4EmStandardPhysics_option4(ver, name), ionScatteringMode_(options.scattering.mode)
{
  SetPhysicsName(name);
  G4CCCurrentChargeConfig::Configure(options);
  G4CCCurrentChargeConfig::Instance().EnsureMessenger();
}

G4EmStandardPhysics_option4_CurrentCharge::G4EmStandardPhysics_option4_CurrentCharge(
  G4int ver, const G4String& name)
  : G4EmStandardPhysics_option4_CurrentCharge(
      ver, IonScatteringMode::Opt4MultipleScattering, name)
{}

G4EmStandardPhysics_option4_CurrentCharge::G4EmStandardPhysics_option4_CurrentCharge(
  G4int ver, IonScatteringMode mode, const G4String& name)
  : G4EmStandardPhysics_option4(ver, name), ionScatteringMode_(mode)
{
  SetPhysicsName(name);
  G4CCCurrentChargeConfig::Instance().EnsureMessenger();
}

void G4EmStandardPhysics_option4_CurrentCharge::ConstructProcess()
{
  G4CCCurrentChargeConfig::Instance().Freeze();
  G4EmStandardPhysics_option4::ConstructProcess();
  G4CCCurrentChargeConfig::Instance().EnsureMessenger();
  ReplaceIonIonisationAndAttachCurrentChargeProcesses();
  if(ionScatteringMode_ == IonScatteringMode::Opt4MultipleScattering) {
    ReplaceNuclearStoppingWithEnergyBudget();
  }
  AttachIonSingleScattering();
}

void G4EmStandardPhysics_option4_CurrentCharge::ReplaceIonIonisationAndAttachCurrentChargeProcesses()
{
  auto* iterator = GetParticleIterator();
  iterator->reset();
  while((*iterator)()) {
    auto* particle = iterator->value();
    if(!G4CCIonChargeState::IsIonDefinition(particle)) { continue; }
    if(G4CCIonChargeState::GetZp(particle) < 29) { continue; }

    auto* pm = particle->GetProcessManager();
    if(pm == nullptr) { continue; }

    RemoveProcessesNamed(pm, "ionIoni");

    if(!HasProcessNamed(pm, "CCIonIoni")) {
      RegisterProcess(new G4CCIonIonisation(), particle);
    }
    if(!HasProcessNamed(pm, "CCIonDeltaElectron")) {
      pm->AddDiscreteProcess(new G4CCIonDeltaElectronProcess());
    }
    if(!HasProcessNamed(pm, "CCIonChargeExchange")) {
      pm->AddDiscreteProcess(new G4CCIonChargeExchangeProcess());
    }
    if(!HasProcessNamed(pm, "CCNonEMChargeReset")) {
      auto* reset = new G4CCNonEMChargeResetProcess();
      pm->AddDiscreteProcess(reset);
      pm->SetProcessOrderingToLast(reset, idxPostStep);
    }
  }

  auto* genericIon = G4GenericIon::GenericIon();
  if(genericIon != nullptr && genericIon->GetProcessManager() != nullptr &&
     !HasProcessNamed(genericIon->GetProcessManager(), "CCIonIoni")) {
    auto* pm = genericIon->GetProcessManager();
    RemoveProcessesNamed(pm, "ionIoni");
    RegisterProcess(new G4CCIonIonisation(), genericIon);
    pm->AddDiscreteProcess(new G4CCIonDeltaElectronProcess());
    pm->AddDiscreteProcess(new G4CCIonChargeExchangeProcess());
    auto* reset = new G4CCNonEMChargeResetProcess();
    pm->AddDiscreteProcess(reset);
    pm->SetProcessOrderingToLast(reset, idxPostStep);
  }
}

void G4EmStandardPhysics_option4_CurrentCharge::AttachIonSingleScattering()
{
  if(ionScatteringMode_ == IonScatteringMode::Opt4MultipleScattering) { return; }
  const G4bool useTables =
    ionScatteringMode_ == IonScatteringMode::TabulatedSingleScattering;
  auto& config = G4CCCurrentChargeConfig::Instance();
  std::shared_ptr<const G4CCElasticTableData> scatteringData;
  if(useTables && config.IsDiracMode()) {
    try {
      scatteringData = G4CCElasticTableData::LoadConfigured(config.Options().scattering.tables);
    } catch(const std::exception& error) {
      G4Exception("G4EmStandardPhysics_option4_CurrentCharge", "CCElasticData001",
                  FatalException, error.what());
      throw;
    }
  }
  auto* genericIon = G4GenericIon::GenericIon();
  std::set<G4ProcessManager*> seen;
  auto* iterator = GetParticleIterator();
  iterator->reset();
  while((*iterator)()) {
    auto* particle = iterator->value();
    if(particle != genericIon &&
       !(particle->GetParticleType() == "nucleus" &&
         particle->GetAtomicNumber() > 2 && !particle->IsShortLived())) { continue; }
    auto* pm = particle->GetProcessManager();
    if(pm == nullptr || !seen.insert(pm).second) { continue; }
    for(G4int i = pm->GetProcessListLength() - 1; i >= 0; --i) {
      auto* process = (*pm->GetProcessList())[i];
      const auto& name = process->GetProcessName();
      if(dynamic_cast<G4VMultipleScattering*>(process) != nullptr ||
         name == "nuclearStopping" || name == "CoulombScat" ||
         name == "CCTabulatedIonScattering") {
        pm->RemoveProcess(i);
      }
    }
    pm->AddDiscreteProcess(new G4CCTabulatedIonScattering(scatteringData, useTables, 0., !config.IsDiracMode()));
  }
}
