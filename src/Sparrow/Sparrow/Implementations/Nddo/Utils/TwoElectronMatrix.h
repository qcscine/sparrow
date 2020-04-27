/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_TWOELECTRONMATRIX_H
#define SPARROW_TWOELECTRONMATRIX_H

#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {

namespace Utils {
enum class ElementType : unsigned;
class DensityMatrix;
class AtomsOrbitalsIndexes;
enum class derivativeType;
} // namespace Utils

namespace Sparrow {

namespace nddo {
class OneCenterIntegralContainer;
class OneCenterTwoElectronIntegrals;
class TwoCenterIntegralContainer;
class ElementParameters;
class AtomicParameters;

namespace multipole {
class Global2c2eMatrix;
}

/*!
 * @brief Class to generate the two-electron matrix G for semi-empirical methods.
 * This class is parallelized with OpenMP.
 */

class TwoElectronMatrix {
 public:
  TwoElectronMatrix(const Utils::ElementTypeCollection& elements, const Utils::DensityMatrix& densityMatrix,
                    const OneCenterIntegralContainer& oneCIntegrals, const TwoCenterIntegralContainer& twoCIntegrals,
                    const ElementParameters& elementPar, const Utils::AtomsOrbitalsIndexes& aoIndexes);
  void initialize();

  void calculate(bool spinPolarized);
  void calculateBlocks();
  void calculateSameAtomBlock(int startIndex, int nAOs, Utils::ElementType el, Eigen::MatrixXd& G,
                              Eigen::MatrixXd& GAlpha, Eigen::MatrixXd& GBeta);
  void calculateDifferentAtomsBlock(int startA, int startB, int nAOsA, int nAOsB, const multipole::Global2c2eMatrix& m,
                                    Eigen::MatrixXd& G, Eigen::MatrixXd& GAlpha, Eigen::MatrixXd& GBeta);
  template<Utils::derivativeType O>
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivativeContainer) const;
  const Eigen::MatrixXd& operator()() const {
    return G_;
  }
  const Eigen::MatrixXd& getMatrix() const {
    return G_;
  }
  const Eigen::MatrixXd& getAlpha() const {
    return GAlpha_;
  }
  const Eigen::MatrixXd& getBeta() const {
    return GBeta_;
  }

 private:
  template<Utils::derivativeType O>
  void addDerivativesForBlock(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivativeContainer, int a, int b,
                              int startA, int startB, int nAOsA, int nAOsB, const multipole::Global2c2eMatrix& m) const;

  bool spinPolarized_;
  const Eigen::MatrixXd &P, &PAlpha_, &PBeta_;
  const OneCenterIntegralContainer& oneCenterIntegrals;
  const TwoCenterIntegralContainer& twoCenterIntegrals;
  const ElementParameters& elementParameters;
  const Utils::AtomsOrbitalsIndexes& aoIndexes_;

  Eigen::MatrixXd G_, GAlpha_, GBeta_;
  const Utils::ElementTypeCollection& elementTypes_;
  int nAOs_;
  int nAtoms_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_TWOELECTRONMATRIX_H
