/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ZeroOrderMatricesCalculator.h"
#include "SKPair.h"
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Utils/Math/DerivOrderEnum.h>
#include <Utils/Typenames.h>
#include <cmath>

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;
using std::sqrt;
namespace {
const double sqrt3 = std::sqrt(3.0);
}

namespace dftb {
ZeroOrderMatricesCalculator::ZeroOrderMatricesCalculator(const Utils::ElementTypeCollection& elements,
                                                         const Utils::PositionCollection& positions,
                                                         const Utils::AtomsOrbitalsIndexes& aoIndexes,
                                                         const DFTBCommon::AtomicParameterContainer& atomicPar,
                                                         const DFTBCommon::DiatomicParameterContainer& diatomicPar,
                                                         const Utils::DensityMatrix& densityMatrix)
  : elements_(elements),
    positions_(positions),
    aoIndexes_(aoIndexes),
    atomicPar_(atomicPar),
    diatomicPar_(diatomicPar),
    densityMatrix_(densityMatrix) {
}

const Utils::MatrixWithDerivatives& ZeroOrderMatricesCalculator::getZeroOrderHamiltonian() const {
  return zeroOrderHamiltonian_;
}

void ZeroOrderMatricesCalculator::calculateOverlap(Utils::derivOrder highestRequiredOrder) {
  constructH0S(highestRequiredOrder);
}

const Utils::MatrixWithDerivatives& ZeroOrderMatricesCalculator::getOverlap() const {
  return overlap_;
}

void ZeroOrderMatricesCalculator::resetOverlap() {
  initializeH0S();
}

void ZeroOrderMatricesCalculator::initializeH0S() {
  // Calculate elements of S and H0 matrices for orbitals on same atom; they don't change at all
  auto nAOs = aoIndexes_.getNAtomicOrbitals();

  Eigen::MatrixXd H0 = Eigen::MatrixXd::Zero(nAOs, nAOs);
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(nAOs, nAOs);
#pragma omp parallel for
  for (int A = 0; A < elements_.size(); A++) {
    int nAOsA = aoIndexes_.getNOrbitals(A);
    int AOindexA = aoIndexes_.getFirstOrbitalIndex(A);
    for (int i = 0; i < nAOsA; i++) {
      for (int j = 0; j < nAOsA; j++) {
        if (i == j) {
          S(AOindexA + i, AOindexA + j) = 1.0;
          H0(AOindexA + i, AOindexA + j) = atomicPar_[Utils::ElementInfo::Z(elements_[A])]->getOrbitalEnergy(i);
        }
        else {
          S(AOindexA + i, AOindexA + j) = 0.0;
          H0(AOindexA + i, AOindexA + j) = 0.0;
        }
      }
    }
  }
  overlap_.setBaseMatrix(S);
  zeroOrderHamiltonian_.setBaseMatrix(H0);
}

void ZeroOrderMatricesCalculator::constructH0S(Utils::derivOrder order) {
  if (order == Utils::derivOrder::zero)
    constructH0S<Utils::derivOrder::zero>();
  else if (order == Utils::derivOrder::one)
    constructH0S<Utils::derivOrder::one>();
  else if (order == Utils::derivOrder::two)
    constructH0S<Utils::derivOrder::two>();
}

template<Utils::derivOrder O>
void ZeroOrderMatricesCalculator::constructH0S() {
  zeroOrderHamiltonian_.setOrder(O);
  overlap_.setOrder(O);

  constructPartOfH0S<O>();

  auto& H0 = zeroOrderHamiltonian_.get<O>();
  auto& S = overlap_.get<O>();
  for (int i = 0; i < aoIndexes_.getNAtomicOrbitals(); i++) {
    for (int j = i + 1; j < aoIndexes_.getNAtomicOrbitals(); j++) {
      S(j, i) = getValueWithOppositeDerivative<O>(S(i, j));
      H0(j, i) = getValueWithOppositeDerivative<O>(H0(i, j));
    }
  }
}

template<Utils::derivOrder O>
void ZeroOrderMatricesCalculator::constructPartOfH0S() {
  using Val = Value3DType<O>;
  auto& H0 = zeroOrderHamiltonian_.get<O>();
  auto& S = overlap_.get<O>();

#pragma omp parallel for // private(me, v, v2, vv, I, val)
  for (int a = 0; a < aoIndexes_.getNAtoms(); a++) {
    Val me[2][9][9]; // Matrix elements; me[0][][] -> overlap; me[1][][] -> hamiltonian
    Val v[3];        // x, y, z
    Val v2[3][3];    // x*x, x*y, y*y, etc.
    Val vv[3];
    Val I[28];
    InterpolationValues<O> val{};

    int nAOsA = aoIndexes_.getNOrbitals(a);
    int AOindexA = aoIndexes_.getFirstOrbitalIndex(a);
    for (int b = a + 1; b < aoIndexes_.getNAtoms(); b++) {
      int nAOsB = aoIndexes_.getNOrbitals(b);
      int AOindexB = aoIndexes_.getFirstOrbitalIndex(b);

      Eigen::Vector3d R = positions_.row(b) - positions_.row(a);
      double dist = R.norm();

      SKPair* parameters;
      if (elements_[a] <= elements_[b])
        parameters = diatomicPar_[Utils::ElementInfo::Z(elements_[a])][Utils::ElementInfo::Z(elements_[b])].get();
      else {
        parameters = diatomicPar_[Utils::ElementInfo::Z(elements_[b])][Utils::ElementInfo::Z(elements_[a])].get();
        R *= -1.0;
      }

      if (parameters->getHS(dist, val) == 0) { // if all values and derivatives are zero
        for (int i = 0; i < nAOsA; i++) {
          for (int j = 0; j < nAOsB; j++) {
            S(AOindexA + i, AOindexB + j) = constant3D<O>(0.0);
            S(AOindexB + j, AOindexA + i) = constant3D<O>(0.0);
            H0(AOindexA + i, AOindexB + j) = constant3D<O>(0.0);
            H0(AOindexB + j, AOindexA + i) = constant3D<O>(0.0);
          }
        }
        continue; // jump to next atom pair
      }
      for (int i = 0; i < parameters->getNIntegrals(); ++i)
        I[i] = get3Dfrom1D<O>(val.derivIntegral[i], R);

      /*
       * NB: notation for integrals
       *
       Ssss=I[0];
       Ssps=I[2];
       Spps=I[4];
       Sppp=I[6];
       Spss=I[8];
       Ssds=I[10];
       Sdss=I[12]; //TODO: Put this after pdp
       Spds=I[14];
       Spdp=I[16];
       Sdps=I[18];
       Sdpp=I[20];
       Sdds=I[22];
       Sddp=I[24];
       Sddd=I[26];

       Hsss=I[1];
       Hsps=I[3];
       Hpps=I[5];
       Hppp=I[7];
       Hpss=I[9];
       Hsds=I[11];
       Hpds=I[15];
       Hpdp=I[17];
       Hdds=I[23];
       Hddp=I[25];
       Hddd=I[27];
       */

      //=======================================================================================================================
      //                                                  s - s
      //=======================================================================================================================
      for (int m = 0; m < 2; m++) {
        me[m][0][0] = I[0 + m];
      }

      if (nAOsA + nAOsB != 2) {
        // Set v[] and v2[]:
        // direction cosines of R and their square
        Val DIM[3] = {toX<O>(R.x()), toY<O>(R.y()), toZ<O>(R.z())};
        auto R2 = toRSquared<O>(R.x(), R.y(), R.z());
        auto RNorm = sqrt(R2);

        for (int i = 0; i < 3; i++) {
          v[i] = DIM[i] / RNorm;
          v2[i][i] = v[i] * v[i];
          for (int j = i - 1; j >= 0; j--)
            v2[i][j] = v2[j][i] = v[i] * v[j];
        }

        for (int m = 0; m < 2; m++) {
          //=======================================================================================================================
          //                                                  s - p
          //=======================================================================================================================

          // s - x,y,z
          for (int j = 0; j < 3; j++)
            me[m][0][j + 1] = I[2 + m] * v[j];

          if (nAOsA != 1 && nAOsB != 1) { // Means there are p-p interactions

            //=======================================================================================================================
            //                                                  p - s
            //=======================================================================================================================

            // Set other s-p interaction: x,y,z - s
            for (int j = 0; j < 3; j++)
              me[m][j + 1][0] = I[8 + m] * v[j];

            //=======================================================================================================================
            //                                                  p - p
            //=======================================================================================================================
            for (int i = 0; i < 3; i++) {
              me[m][i + 1][i + 1] = I[4 + m] * v2[i][i] + I[6 + m] * (constant3D<O>(1.0) - v2[i][i]);
              for (int j = i + 1; j < 3; j++)
                me[m][j + 1][i + 1] = me[m][i + 1][j + 1] = v2[i][j] * (I[4 + m] - I[6 + m]);
            }
          }
        }

        if (nAOsA == 9 || nAOsB == 9) { // there are d orbitals
          // Expressions that are often needed and their derivative vectors
          auto alpha = v2[0][0] + v2[1][1];
          auto beta = v2[0][0] - v2[1][1];
          auto z2halpha = v2[2][2] - 0.5 * alpha;
          vv[0] = v2[0][1]; // xy
          vv[1] = v2[1][2]; // yz
          vv[2] = v2[0][2]; // xz

          for (int m = 0; m < 2; m++) {
            //=======================================================================================================================
            //                                                  s - d
            //=======================================================================================================================

            // dxy, dyz, dxz with s
            for (int i = 0; i < 3; i++)
              me[m][0][i + 4] = I[10 + m] * sqrt3 * vv[i];

            // dx2-y2 with s
            me[m][0][7] = I[10 + m] * (0.5 * sqrt3) * beta;

            // d2z2-r2 with s
            me[m][0][8] = I[10 + m] * z2halpha;
          }

          if (nAOsA + nAOsB >= 13) { // i.e. there are d and p orbitals
            auto xyz = v2[0][1] * v[2];

            for (int m = 0; m < 2; m++) {
              //=======================================================================================================================
              //                                                  p - d
              //=======================================================================================================================

              // complementary: x-yz, y-xz, z-xy
              for (int i = 0; i < 3; i++)
                me[m][i + 1][(i + 1) % 3 + 4] = xyz * (sqrt3 * I[14 + m] - 2 * I[16 + m]);

              // x-xy, x-xz, y-xy, y-yz, z-xz, z-yz
              for (int i = 0; i < 3; i++) {
                for (int k = 0; k < 2; k++) {
                  int j = (i - k + 3) % 3;
                  int l = (4 - i - j) % 3;
                  me[m][i + 1][j + 4] = sqrt3 * v[i] * vv[j] * I[14 + m] + (v[l] - 2 * v[i] * vv[j]) * I[16 + m];
                }
              }

              // x,y,z - dx2-y2
              for (int i = 0; i < 3; i++)
                me[m][i + 1][7] = v[i] * beta * (0.5 * sqrt3 * I[14 + m] - I[16 + m]); // not complete yet, see below
              // Add missing terms
              me[m][1][7] += v[0] * I[16 + m];
              me[m][2][7] -= v[1] * I[16 + m];

              // x,y - d3z2-r2
              for (int i = 0; i < 2; i++)
                me[m][i + 1][8] = v[i] * z2halpha * I[14 + m] - sqrt3 * v[i] * v2[2][2] * I[16 + m];

              // z - d3z2-r2
              me[m][3][8] = v[2] * z2halpha * I[14 + m] + sqrt3 * v[2] * alpha * I[16 + m];
            }

            if (nAOsA + nAOsB == 18) { // Two d orbitals

              for (int m = 0; m < 2; m++) {
                //=======================================================================================================================
                //                                                  d - s
                //=======================================================================================================================

                // dxy, dyz, dxz with s
                for (int i = 0; i < 3; i++)
                  me[m][i + 4][0] = I[12 + m] * sqrt3 * vv[i];

                // dx2-y2 with s
                me[m][7][0] = I[12 + m] * 0.5 * sqrt3 * beta;

                // d2z2-r2 with s
                me[m][8][0] = I[12 + m] * z2halpha;

                //=======================================================================================================================
                //                                                  d - p
                //=======================================================================================================================

                // complementary: x-yz, y-xz, z-xy
                for (int i = 0; i < 3; i++)
                  me[m][(i + 1) % 3 + 4][i + 1] = xyz * (sqrt3 * I[18 + m] - 2 * I[20 + m]);

                // x-xy, x-xz, y-xy, y-yz, z-xz, z-yz
                for (int i = 0; i < 3; i++) {
                  for (int k = 0; k < 2; k++) {
                    int j = (i - k + 3) % 3;
                    int l = (4 - i - j) % 3;
                    me[m][j + 4][i + 1] = sqrt3 * v[i] * vv[j] * I[18 + m] + (v[l] - 2 * v[i] * vv[j]) * I[20 + m];
                  }
                }

                // x,y,z - dx2-y2
                for (int i = 0; i < 3; i++)
                  me[m][7][i + 1] = v[i] * beta * (0.5 * sqrt3 * I[18 + m] - I[20 + m]); // not complete yet, see below
                // Add missing terms
                me[m][7][1] += v[0] * I[20 + m];
                me[m][7][2] -= v[1] * I[20 + m];

                // x,y - d3z2-r2
                for (int i = 0; i < 2; i++)
                  me[m][8][i + 1] = v[i] * z2halpha * I[18 + m] - sqrt3 * v[i] * v2[2][2] * I[20 + m];

                // z - d3z2-r2
                me[m][8][3] = v[2] * z2halpha * I[18 + m] + sqrt3 * v[2] * alpha * I[20 + m];

                //=======================================================================================================================
                //                                                  d - d
                //=======================================================================================================================

                for (int i = 0; i < 3; i++) {
                  // xy-xy, yz-yz, xz-xz
                  int i1 = i, i2 = (i + 1) % 3, i3 = (i + 2) % 3;
                  me[m][i + 4][i + 4] = vv[i] * vv[i] * (3 * I[22 + m] - 4 * I[24 + m] + I[26 + m]) +
                                        (v2[i1][i1] + v2[i2][i2]) * I[24 + m] + v2[i3][i3] * I[26 + m];
                  // xy-xz, xy-yz, yz-xz
                  for (int j = i + 1; j < 3; j++)
                    me[m][i + 4][j + 4] = me[m][j + 4][i + 4] = vv[i] * vv[j] * (3 * I[22 + m] - 4 * I[24 + m] + I[26 + m]) +
                                                                vv[3 - i - j] * (I[24 + m] - I[26 + m]);
                }

                for (int i = 0; i < 3; i++) {
                  double factor = (i == 0 ? 0.0 : (i == 1 ? 1.0 : -1.0));
                  // xy - x2-y2, xz - x2-y2, yz - x2-y2
                  me[m][4 + i][7] = me[m][7][4 + i] = vv[i] * beta * (1.5 * I[22 + m] - 2 * I[24 + m] + 0.5 * I[26 + m]) +
                                                      factor * vv[i] * (-I[24 + m] + I[26 + m]);
                }

                // xy - 3z2-r2
                me[m][4][8] = me[m][8][4] = sqrt3 * (vv[1] * vv[2] * (I[22 + m] - 2 * I[24 + m] + 0.5 * I[26 + m]) -
                                                     0.5 * vv[0] * alpha * I[22 + m] + 0.5 * vv[0] * I[26 + m]);

                for (int i = 1; i < 3; i++) {
                  // yz - 3z2-r2, xz - 3z2-r2
                  me[m][4 + i][8] = me[m][8][4 + i] =
                      sqrt3 * (vv[i] * (alpha * (-0.5 * I[22 + m] + I[24 + m] - 0.5 * I[26 + m]) +
                                        v2[2][2] * (I[22 + m] - I[24 + m])));
                }

                // x2-y2 - x2-y2
                me[m][7][7] = beta * beta * (0.75 * I[22 + m] - I[24 + m] + 0.25 * I[26 + m]) + alpha * I[24 + m] +
                              v2[2][2] * I[26 + m];

                // x2-y2 - 3z2-r2
                me[m][7][8] = me[m][8][7] = sqrt3 * beta *
                                            (v2[2][2] * (0.5 * I[22 + m] - I[24 + m] + 0.25 * I[26 + m]) -
                                             0.25 * alpha * I[22 + m] + 0.25 * I[26 + m]);

                // 3z2-r2 - 3z2-r2
                me[m][8][8] = z2halpha * z2halpha * I[22 + m] + 3 * v[2] * v[2] * alpha * I[24 + m] +
                              0.75 * alpha * alpha * I[26 + m];
              }

            } // End d-d
          }   // End d-p

        } // End d
      }

      // Copy arrays into S and H matrices
      if (elements_[a] <= elements_[b]) {
        for (int i = 0; i < nAOsA; i++) {
          for (int j = 0; j < nAOsB; j++) {
            S(AOindexA + i, AOindexB + j) = me[0][i][j];
            H0(AOindexA + i, AOindexB + j) = me[1][i][j];
          }
        }
      }
      else {
        for (int i = 0; i < nAOsA; i++) {
          for (int j = 0; j < nAOsB; j++) {
            S(AOindexA + i, AOindexB + j) = getValueWithOppositeDerivative<O>(me[0][j][i]);
            H0(AOindexA + i, AOindexB + j) = getValueWithOppositeDerivative<O>(me[1][j][i]);
          }
        }
      }
    }
  }
}

void ZeroOrderMatricesCalculator::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives,
    const Eigen::MatrixXd& overlapDerivativeMultiplier) const {
  addDerivativesImpl<Utils::derivativeType::first>(derivatives, overlapDerivativeMultiplier);
}

void ZeroOrderMatricesCalculator::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives,
    const Eigen::MatrixXd& overlapDerivativeMultiplier) const {
  addDerivativesImpl<Utils::derivativeType::second_atomic>(derivatives, overlapDerivativeMultiplier);
}

void ZeroOrderMatricesCalculator::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives,
    const Eigen::MatrixXd& overlapDerivativeMultiplier) const {
  addDerivativesImpl<Utils::derivativeType::second_full>(derivatives, overlapDerivativeMultiplier);
}

template<Utils::derivativeType O>
void ZeroOrderMatricesCalculator::addDerivativesImpl(DerivativeContainerType<O>& derivatives,
                                                     const Eigen::MatrixXd& overlapDerivativeMultiplier) const {
  auto nAtoms = elements_.size();

  Value3DType<UnderlyingOrder<O>> der;
  DerivativeType<O> derivative;
  derivative.setZero();
#pragma omp parallel for firstprivate(derivative) private(der)
  for (int a = 0; a < nAtoms; ++a) {
    int nAOsA = aoIndexes_.getNOrbitals(a);
    int AOindexA = aoIndexes_.getFirstOrbitalIndex(a);

    for (int b = a + 1; b < nAtoms; b++) {
      int nAOsB = aoIndexes_.getNOrbitals(b);
      int AOindexB = aoIndexes_.getFirstOrbitalIndex(b);

      der = constant3D<UnderlyingOrder<O>>(0);
      for (int i = 0; i < nAOsA; i++) {
        for (int j = 0; j < nAOsB; j++) {
          double Pel = densityMatrix_.restricted(AOindexA + i, AOindexB + j);
          double Wel = overlapDerivativeMultiplier(AOindexA + i, AOindexB + j);
          der += 2 * (Pel * zeroOrderHamiltonian_.get<UnderlyingOrder<O>>()(AOindexA + i, AOindexB + j) -
                      Wel * overlap_.get<UnderlyingOrder<O>>()(AOindexA + i, AOindexB + j));
        }
      }
      derivative = getDerivativeFromValueWithDerivatives<O>(der);
#pragma omp critical
      { addDerivativeToContainer<O>(derivatives, a, b, derivative); }
    }
  }
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
