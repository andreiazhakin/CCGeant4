#ifndef TransportPhysicsList_hh
#define TransportPhysicsList_hh

#include "G4VModularPhysicsList.hh"

#include <string>

class PhysicsList : public G4VModularPhysicsList
{
public:
  PhysicsList(const std::string& mode, const std::string& scattering,
              const std::string& tables, const std::string& elasticTables,
              double cut);

  void SetCuts() override;
};

#endif
