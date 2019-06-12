/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_SKPAIR_H
#define SPARROW_DFTB_SKPAIR_H

#include "RepulsionParameters.h"
#include "Utils/Math/DerivOrderEnum.h"
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace dftb {
class SKAtom;

template<Utils::derivOrder O>
struct InterpolationValues {
  Utils::AutomaticDifferentiation::Value1DType<O> derivIntegral[28];
  Utils::AutomaticDifferentiation::Value1DType<O> C[28][8];
  Utils::AutomaticDifferentiation::Value1DType<O> D[28][8];
  double xa[8];
  double ya[28][8];
};

// Class containing the DFTB0 parameters of a pair of atoms
class SKPair {
 public:
  SKPair(const std::string& atom1, const std::string& atom2, SKAtom* atomicParameters1, SKAtom* atomicParameters2,
         const std::string& path);
  void complete(SKPair* p);

  template<Utils::derivOrder O>
  int getHS(double dist, InterpolationValues<O>& val) const;
  template<Utils::derivOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> getRepulsion(double const& r) const;
  void getGammaTerms(double& t1a, double& t1b, double& t2a, double& t2b);
  void getGammaDerTerms(double& dtab1a, double& dtab1b, double& dtab2a, double& dtab2b, double& dtba1a, double& dtba1b,
                        double& dtba2a, double& dtba2b, double& dtadr, double& dtbdr);
  void precalculateGammaTerms();
  int getNIntegrals() const {
    return nIntegrals;
  }
  const dftb::RepulsionParameters& getRepulsionParameters() const {
    return repulsion_;
  }

 private:
  // Private methods
  template<Utils::derivOrder O>
  int getHSIntegral(InterpolationValues<O>& val, double dist) const;
  template<Utils::derivOrder O>
  void interpolate(InterpolationValues<O>& val, double x, int start) const;
  void precompute5Extrapolation();
  void readSKFfile(const std::string& atom1, const std::string& atom2, SKAtom* atomicParameters1,
                   SKAtom* AtomParameters2, const std::string& path);

  // Private member variables
  SKAtom* atomType1;
  SKAtom* atomType2;
  double gridDist;
  int nGridPoints;
  double rMax;                        // Distance after which integral is zero
  std::vector<std::vector<double>> M; // Matrix containing the integrals

  dftb::RepulsionParameters repulsion_;

  double distFudge; // Distance from last integral value to zero

  // Members used for the interpolation
  std::vector<int> integralIndexes;
  std::vector<double> extrC3, extrC4, extrC5; // precomputed coefficients for integral extrapolation
  int nInter, nInterRight, nInterLeft;        // Number of points to consider for interpolation
  int nIntegrals;                             // Number of integrals that need be calculated for the pair
  double deltaR;                              // Used to calculate numerical derivatives

  // Members for efficient gamma calculation
  double g1a, g1b, g2a, g2b;
  double dgab1a, dgab1b, dgab2a, dgab2b, dgba1a, dgba1b, dgba2a, dgba2b, dgadr, dgbdr;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif
