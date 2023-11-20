/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "PM6PairwiseRepulsion.h"

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace nddo {

PM6PairwiseRepulsion::PM6PairwiseRepulsion(const AtomicParameters& A, const AtomicParameters& B, const PM6DiatomicParameters& AB)
  : cOfAdditiveTerm(0.0000207389327971913490822259), // NB: The value of c = 1e-8 described in the paper is transformed
                                                     // from Angstrom^12 to bohr^12.
    exponentCoefficient(0.00001244878365801758178693335), // NB: The value of 0.0003 described in the paper is
                                                          // transformed from Angstrom^-5 to bohr^-5.
    furtherExponentCC(3.1644797213016), // NB: The value of 5.98 described in the paper is transformed from Angstrom^-1
                                        // to bohr^-1.
    distanceSiO(5.480205761238679760340424), // NB: The value of 2.9 described in the paper is transformed from Angstrom
                                             // to bohr.
    factorCC(9.28),
    factorSiO(-0.0007),
    pA_(A),
    pB_(B),
    pAB_(AB) {
}

void PM6PairwiseRepulsion::calculate(const Eigen::Vector3d& R, Utils::DerivativeOrder order) {
  double Rabs = R.norm();
  if (order == Utils::DerivativeOrder::Zero) {
    repulsionEnergy_ = calculateRepulsion<Utils::DerivativeOrder::Zero>(Rabs);
  }
  else if (order == Utils::DerivativeOrder::One) {
    First1D rep = calculateRepulsion<Utils::DerivativeOrder::One>(Rabs);
    repulsionEnergy_ = rep.value();
    repulsionGradient_ = Utils::Gradient(get3Dfrom1D<Utils::DerivativeOrder::One>(rep, R).derivatives());
  }
  else if (order == Utils::DerivativeOrder::Two) {
    Second1D rep = calculateRepulsion<Utils::DerivativeOrder::Two>(Rabs);
    repulsionEnergy_ = rep.value();
    repulsionHessian_ = get3Dfrom1D<Utils::DerivativeOrder::Two>(rep, R);
  }
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
