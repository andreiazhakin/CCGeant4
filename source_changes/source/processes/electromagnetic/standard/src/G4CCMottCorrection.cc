#include "G4CCMottCorrection.hh"

#include "G4MottData.hh"

#include <algorithm>
#include <cmath>

namespace
{
  inline G4double Clamp01(G4double x)
  {
    if(x < 0.0) { return 0.0; }
    if(x > 1.0) { return 1.0; }
    return x;
  }

  inline G4double EvaluateMottPolynomial(const G4double coeff[5],
                                        G4double deltaEnergy, G4double tMax)
  {
    if(deltaEnergy <= 0.0 || tMax <= 0.0) { return 0.0; }
    const G4double x = Clamp01(deltaEnergy/tMax);
    const G4double y = std::sqrt(2.0*x);

    G4double yPow = 1.0;
    G4double result = 0.0;
    for(G4int j = 0; j <= 4; ++j) {
      result += coeff[j]*yPow;
      yPow *= y;
    }
    return std::max(0.0, result);
  }
}

G4bool G4CCMottCorrection::IsAvailable(G4int zProjectile)
{
  return zProjectile >= 1 && zProjectile <= 92;
}

G4bool G4CCMottCorrection::PolynomialCoefficients(G4int zProjectile,
                                                         G4double beta,
                                                         G4double coeff[5])
{
  if(coeff == nullptr) { return false; }
  for(G4int j = 0; j <= 4; ++j) { coeff[j] = 0.0; }

  if(!IsAvailable(zProjectile) || beta <= 0.0) { return false; }

  const G4double db = beta - betaBar;
  G4double betaPow[6];
  betaPow[0] = 1.0;
  for(G4int k = 1; k < 6; ++k) { betaPow[k] = betaPow[k-1]*db; }

  for(G4int j = 0; j <= 4; ++j) {
    G4double aj = 0.0;
    for(G4int k = 0; k <= 5; ++k) {
      aj += fMottCoef[zProjectile][j][k]*betaPow[k];
    }
    coeff[j] = aj;
  }
  return true;
}

G4double G4CCMottCorrection::FullStoppingCorrection(G4int zProjectile,
                                                    G4double beta)
{
  if(!std::isfinite(beta) || beta <= 0.0 || beta >= 1.0) { return 0.0; }

  G4double coeff[5];
  if(!PolynomialCoefficients(zProjectile, beta, coeff)) { return 0.0; }
  for(const G4double value : coeff) {
    if(!std::isfinite(value)) { return 0.0; }
  }

  const G4double sqrt2 = std::sqrt(2.0);
  const G4double result = beta*beta + 2.0*sqrt2*coeff[1] + 2.0*coeff[2]
    + (4.0*sqrt2/3.0)*coeff[3] + 2.0*coeff[4];
  return std::isfinite(result) ? result : 0.0;
}

G4double G4CCMottCorrection::RMott(G4int zProjectile, G4double beta,
                                   G4double deltaEnergy, G4double tMax)
{
  if(deltaEnergy <= 0.0 || tMax <= 0.0) { return 0.0; }

  G4double coeff[5];
  if(!PolynomialCoefficients(zProjectile, beta, coeff)) { return 0.0; }

  return EvaluateMottPolynomial(coeff,deltaEnergy,tMax);
}

G4double G4CCMottCorrection::IntegratedRMottOverT2(G4int zProjectile,
                                                   G4double beta,
                                                   G4double tMin,
                                                   G4double tMax)
{
  if(tMax <= tMin || tMin <= 0.0 || tMax <= 0.0) { return 0.0; }

  G4double coeff[5];
  if(!PolynomialCoefficients(zProjectile, beta, coeff)) { return 0.0; }

  G4double result = 0.0;
  for(G4int j = 0; j <= 4; ++j) {
    const G4double power = 0.5*static_cast<G4double>(j);
    const G4double scale = std::pow(2.0/tMax, power);

    G4double primitiveDifference = 0.0;
    if(j == 2) {
      primitiveDifference = std::log(tMax/tMin);
    } else {
      const G4double exponent = power - 1.0;
      primitiveDifference =
        (std::pow(tMax, exponent) - std::pow(tMin, exponent))/exponent;
    }
    result += coeff[j]*scale*primitiveDifference;
  }

  return std::max(0.0, result);
}

G4double G4CCMottCorrection::G4BornShape(G4double beta2, G4double deltaEnergy,
                                         G4double tMax, G4double totalEnergy2,
                                         G4double spin)
{
  if(deltaEnergy <= 0.0 || tMax <= 0.0) { return 0.0; }
  G4double value = 1.0 - beta2*deltaEnergy/tMax;
  if(spin == 0.5 && totalEnergy2 > 0.0) {
    value += 0.5*deltaEnergy*deltaEnergy/totalEnergy2;
  }
  return std::max(0.0, value);
}

G4double G4CCMottCorrection::CMottVsG4(G4int zProjectile, G4double beta,
                                       G4double deltaEnergy, G4double tMax,
                                       G4double beta2, G4double totalEnergy2,
                                       G4double spin)
{
  const G4double born = G4BornShape(beta2, deltaEnergy, tMax, totalEnergy2, spin);
  if(born <= 0.0) { return 1.0; }
  return RMott(zProjectile, beta, deltaEnergy, tMax)/born;
}

G4double G4CCMottCorrection::MaxRMott(G4int zProjectile, G4double beta,
                                      G4double tMin, G4double tMax)
{
  if(!IsAvailable(zProjectile) || beta <= 0.0 || tMax <= tMin || tMin <= 0.0) {
    return 1.0;
  }

  G4double coeff[5];
  if(!PolynomialCoefficients(zProjectile, beta, coeff)) { return 1.0; }

  G4double maxValue = 1.0;
  constexpr G4int n = 512;
  const G4double u0 = std::log(tMin);
  const G4double u1 = std::log(tMax);
  for(G4int i = 0; i <= n; ++i) {
    const G4double u = u0 + (u1 - u0)*static_cast<G4double>(i)/static_cast<G4double>(n);
    const G4double t = std::exp(u);
    maxValue = std::max(maxValue, EvaluateMottPolynomial(coeff,t,tMax));
  }

  return 1.05*maxValue;
}
