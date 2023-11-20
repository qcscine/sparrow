/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_SKPAIR_H
#define SPARROW_DFTB_SKPAIR_H

#include "RepulsionParameters.h"
#include "Sparrow/Implementations/Dftb/Utils/SkfParser.h"
#include "Utils/Math/DerivOrderEnum.h"
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace dftb {
class SKAtom;

template<Utils::DerivativeOrder O>
struct InterpolationValues {
  std::array<Utils::AutomaticDifferentiation::Value1DType<O>, 28> derivIntegral;
  std::array<std::array<Utils::AutomaticDifferentiation::Value1DType<O>, 8>, 28> C;
  std::array<std::array<Utils::AutomaticDifferentiation::Value1DType<O>, 8>, 28> D;
  std::array<double, 8> xa;
  std::array<std::array<double, 8>, 28> ya;
};

// Class containing the DFTB0 parameters of a pair of atoms
class SKPair {
 public:
  SKPair(SKAtom* atomicParameters1, SKAtom* atomicParameters2, SkfData data);

  void complete(SKPair* p);

  template<Utils::DerivativeOrder O>
  int getHS(double dist, InterpolationValues<O>& val) const;
  template<Utils::DerivativeOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> getRepulsion(double const& r) const;

  struct GammaTerms {
    double g1a, g1b, g2a, g2b;
  };
  const GammaTerms& getGammaTerms() const;

  struct GammaDerivativeTerms {
    double dgab1a, dgab1b, dgab2a, dgab2b, dgba1a, dgba1b, dgba2a, dgba2b, dgadr, dgbdr;
  };
  const GammaDerivativeTerms& getGammaDerTerms() const;

  void precalculateGammaTerms();
  int getNIntegrals() const {
    return nIntegrals;
  }
  const dftb::RepulsionParameters& getRepulsionParameters() const {
    return repulsion_;
  }

 private:
  // Private methods
  template<Utils::DerivativeOrder O>
  int getHSIntegral(InterpolationValues<O>& val, double dist) const;
  template<Utils::DerivativeOrder O>
  void interpolate(InterpolationValues<O>& val, double x, int start) const;
  void precompute5Extrapolation();

  //! Distance from last integral value to zero
  static constexpr double distFudge = 1.0;

  //! Number of points to consider for interpolation
  static constexpr int nInter = 8;
  static constexpr int nInterRight = 4;
  static constexpr int nInterLeft = 4;

  //! Step length in num. derivative finite difference calc. in interpolation
  static constexpr double deltaR = 1e-5;

  //! Order of integrals as we need it (s orbitals first, then p, etc.)
  static constexpr std::array<int, 28> integralIndexes = {
      {19, 9, 18, 8, 15, 5, 16, 6, 27, 23, 17, 7, 26, 22, 13, 3, 14, 4, 24, 20, 25, 21, 10, 0, 11, 1, 12, 2}};

  // Private member variables
  SKAtom* atomType1;
  SKAtom* atomType2;
  double gridDist;
  int nGridPoints;
  double rMax;                                       // Distance after which integral is zero
  std::array<std::vector<double>, 28> integralTable; // Matrix containing the integrals

  dftb::RepulsionParameters repulsion_;

  // Members used for the interpolation
  std::vector<double> extrC3, extrC4, extrC5; // precomputed coefficients for integral extrapolation
  int nIntegrals;                             // Number of integrals that need be calculated for the pair

  // Members for efficient gamma calculation
  GammaTerms gamma;
  GammaDerivativeTerms gammaDerivative;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif
