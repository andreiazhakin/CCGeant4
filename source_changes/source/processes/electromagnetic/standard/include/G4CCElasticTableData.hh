#pragma once
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace tabulated {
struct EnergyBracket {
  bool valid = false;
  std::size_t lo = 0, hi = 0;
  double w = 0., sigma_mb = 0.;
};

struct NodeView {
  std::size_t size = 0;
  double energy_MeVu = 0., sigma_mb = 0., As = 0., C_mb = 0., p_cm_MeV = 0.;
  const double *x = nullptr, *R = nullptr, *CDF = nullptr,
               *survival = nullptr, *cell_mb = nullptr;
};

class ProjectileTable final {
 public:
  int projectile_Z = 0, projectile_A = 0, target_Z = 0, target_A = 0;
  double target_mass_MeV = 0., projectile_mass_MeV = 0.;
  std::size_t EnergyCount() const noexcept { return energy_count_; }
  NodeView Node(std::size_t i) const;
  EnergyBracket Lookup(double energy_MeVu) const;
  std::size_t SelectNode(const EnergyBracket&, double uniform) const;
  double SampleNodeX(std::size_t node, double uniform) const;
  double SampleX(const EnergyBracket&, double uniform_node, double uniform_x) const;
  double InverseNode(std::size_t node, double probability, bool survival = false) const;
  double CumulativeNode(std::size_t node, double x, bool survival = false) const;
  double Cumulative(const EnergyBracket&, double x, bool survival = false) const;

 private:
  friend class G4CCElasticTableData;
  std::size_t energy_count_ = 0, point_count_ = 0;
  const double *energies_ = nullptr, *sigma_ = nullptr, *As_ = nullptr,
               *p_cm_ = nullptr, *C_ = nullptr, *x_ = nullptr, *R_ = nullptr,
               *CDF_ = nullptr, *survival_ = nullptr, *cell_ = nullptr;
  const std::int64_t* offsets_ = nullptr;
};

class G4CCElasticTableData final {
 public:
  static std::shared_ptr<const G4CCElasticTableData> Load(
      const std::string& directory,
      const std::vector<std::pair<int, int>>& projectiles = {});
  static std::shared_ptr<const G4CCElasticTableData> LoadFromEnvironment();
  static std::shared_ptr<const G4CCElasticTableData> LoadConfigured(const std::string& path);
  ~G4CCElasticTableData();
  const ProjectileTable* FindByZ(int projectileZ, int targetZ) const;
  const ProjectileTable* Find(int zp, int ap, int zt, int at) const;
  std::size_t LoadedCount() const noexcept;
  std::size_t IndexedCount() const noexcept;
  std::uint64_t MappedBytes() const noexcept;
  G4CCElasticTableData(const G4CCElasticTableData&) = delete;
  G4CCElasticTableData& operator=(const G4CCElasticTableData&) = delete;
 private:
  G4CCElasticTableData();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
using ElasticTableData = G4CCElasticTableData;
}

using G4CCElasticTableData = tabulated::G4CCElasticTableData;
