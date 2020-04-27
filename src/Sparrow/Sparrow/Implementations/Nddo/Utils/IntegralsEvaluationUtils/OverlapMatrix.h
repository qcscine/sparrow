/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_OVERLAPMATRIX_H
#define SPARROW_OVERLAPMATRIX_H

#include "AtomPairOverlap.h"
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/Scf/MethodInterfaces/OverlapCalculator.h>
#include <Utils/Typenames.h>

namespace Scine {

namespace Utils {
class AtomsOrbitalsIndexes;
} // namespace Utils

namespace Sparrow {
namespace nddo {
class ElementParameters;

/**
 * @brief This class computes the whole overlap matrix and returns it in *lower* diagonal form.
 * The basis function overlap, as well as its first and second order derivatives with respect to the nuclear cartesian
 * coordinates is calculated. It inherits from OverlapCalculator in order to make this class compatible with its
 * polymorphic useage.
 */

class OverlapMatrix : public Utils::OverlapCalculator {
 public:
  //! @brief Constructor
  OverlapMatrix(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                const Utils::AtomsOrbitalsIndexes& aoIndexes, const ElementParameters& elementParameters);
  //! @brief Function calculating the overlap between the AO basis functions up to the desired derivative order
  void calculateOverlap(Utils::derivOrder highestRequiredOrder) override;
  //! @brief Getter for the overlap matrix with its derivatives.
  const Utils::MatrixWithDerivatives& getOverlap() const override;
  //! @brief (Re)-initializes the overlap matrix with its derivatives.
  void reset() override;

 private:
  Utils::MatrixWithDerivatives S_;
  const Utils::ElementTypeCollection& elementTypes_;
  const Utils::PositionCollection& positions_;
  const Utils::AtomsOrbitalsIndexes& aoIndexes_;
  const ElementParameters& elementParameters_;
  AtomPairOverlap<Utils::derivOrder::one> pairOverlapFirstOrder_;
  AtomPairOverlap<Utils::derivOrder::zero> pairOverlapZeroOrder_;
  AtomPairOverlap<Utils::derivOrder::two> pairOverlapSecondOrder_;
  int nAOs_ = 0;
  int nAtoms_ = 0;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_OVERLAPMATRIX_H
