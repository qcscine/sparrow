/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ONEELECTRONMATRIX_H
#define SPARROW_ONEELECTRONMATRIX_H

#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {

namespace Utils {
class MatrixWithDerivatives;
class AtomsOrbitalsIndexes;
} // namespace Utils

namespace Sparrow {

namespace nddo {
class AtomicParameters;
class ElementParameters;
class TwoCenterIntegralContainer;

/**
 * @brief This class generates the one-electron matrix H for semi-empirical methods.
 */
class OneElectronMatrix {
 public:
  //! @brief Constructor.
  OneElectronMatrix(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                    const Eigen::MatrixXd& densityMatrix, const TwoCenterIntegralContainer& twoCIntegrals,
                    const ElementParameters& elementPar, const Utils::AtomsOrbitalsIndexes& aoIndexes);
  //! @brief Initializes one-electron matrix H and the used data structures.
  void initialize();
  //! @brief Fills the one-electron matrix H by reading and calculating the required integrals.
  void calculate(const Utils::MatrixWithDerivatives& S);
  //! @brief Calculates all the blocks on the same atoms.
  void calculateSameAtomBlocks();
  //! @brief Calculates a specific block on an atom.
  void calculateSameAtomBlock(int a, int startIndex, int nAOs);
  //! @brief Calculates all the blocks on different atom pairs.
  void calculateDifferentAtomsBlocks(const Utils::MatrixWithDerivatives& S);
  //! @brief Calculates a block between a specific atom pair.
  void calculateDifferentAtomsBlock(int startRow, int startCol, const AtomicParameters& pA, const AtomicParameters& pB,
                                    const Utils::MatrixWithDerivatives& S);
  //! @brief Calculates the derivative contribution up to the order \tparam O.
  template<Utils::derivativeType O>
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivativeContainer,
                      const Utils::MatrixWithDerivatives& S) const;
  //! @brief Getter for the one-electron matrix H.
  const Eigen::MatrixXd& operator()() const {
    return H_;
  }
  const Eigen::MatrixXd& getMatrix() const {
    return H_;
  }

 private:
  template<Utils::derivativeType O>
  void addDerivativesContribution1(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivativeContainer,
                                   int a, int startIndex, int nAOs) const;
  template<Utils::derivativeType O>
  void addDerivativesContribution2(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivativeContainer,
                                   int a, int b, int indexA, int indexB, int nAOsA, int nAOsB,
                                   const Utils::MatrixWithDerivatives& S) const;

  const Eigen::MatrixXd& P;
  const TwoCenterIntegralContainer& twoCenterIntegrals;
  const ElementParameters& elementParameters;
  const Utils::AtomsOrbitalsIndexes& aoIndexes_;

  int nAOs_ = 0;
  int nAtoms_ = 0;
  Eigen::MatrixXd H_;
  const Utils::ElementTypeCollection& elementTypes_;
  const Utils::PositionCollection& positions_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ONEELECTRONMATRIX_H
