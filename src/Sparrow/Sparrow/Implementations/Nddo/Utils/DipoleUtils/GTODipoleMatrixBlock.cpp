/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#define _USE_MATH_DEFINES
#include "GTODipoleMatrixBlock.h"
#include "AnalyticalDipoleIntegralOverGTOsCalculator.h"
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/AnalyticalDipoleIntegralOverGTOsCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/GeneralTypes.h>
#include <Utils/DataStructures/GtoExpansion.h>
#include <chrono>
#include <cmath>

namespace Scine {
namespace Sparrow {

GTODipoleMatrixBlock::GTODipoleMatrixBlock() {
  PminusA_ = {};
  PminusB_ = {};
  PminusC_ = {};
  overlapIntegral_ = {};
  dipoleIntegral_ = {};
  for (auto& dCB : dipoleComponentsBlocks_)
    dCB = Eigen::Matrix<double, 6, 6>::Zero(6, 6);

  method_ = IntegralMethod::ObaraSaika;

  AOMomenta_[0] = nddo::AngularMomentum(0, 0, 0);
  AOMomenta_[1] = nddo::AngularMomentum(1, 0, 0);
  AOMomenta_[2] = nddo::AngularMomentum(0, 1, 0);
  AOMomenta_[3] = nddo::AngularMomentum(0, 0, 1);
  AOMomenta_[4] = nddo::AngularMomentum(1, 0, 1);
  AOMomenta_[5] = nddo::AngularMomentum(0, 1, 1);
  AOMomenta_[6] = nddo::AngularMomentum(1, 1, 0);
  AOMomenta_[7] = nddo::AngularMomentum(0, 0, 2); // z2
  AOMomenta_[8] = nddo::AngularMomentum(0, 2, 0); // y2
  AOMomenta_[9] = nddo::AngularMomentum(2, 0, 0); // x2

  AOIndexes_[0] = static_cast<int>(nddo::GeneralTypes::orb_t::s);
  AOIndexes_[1] = static_cast<int>(nddo::GeneralTypes::orb_t::x) - 1;
  AOIndexes_[2] = static_cast<int>(nddo::GeneralTypes::orb_t::y) - 1;
  AOIndexes_[3] = static_cast<int>(nddo::GeneralTypes::orb_t::z) - 1;
  AOIndexes_[4] = static_cast<int>(nddo::GeneralTypes::orb_t::xz) - 4;
  AOIndexes_[5] = static_cast<int>(nddo::GeneralTypes::orb_t::yz) - 4;
  AOIndexes_[6] = static_cast<int>(nddo::GeneralTypes::orb_t::xy) - 4;
  AOIndexes_[7] = static_cast<int>(nddo::GeneralTypes::orb_t::z2) - 4;
  AOIndexes_[8] = static_cast<int>(nddo::GeneralTypes::orb_t::x2y2) - 4;
  AOIndexes_[9] = 0;
}

void GTODipoleMatrixBlock::calculateSingleGTFPair(int GaussTypeFunctionIndexA, int GaussTypeFunctionIndexB,
                                                  const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB,
                                                  const Eigen::RowVector3d& Ra, const Eigen::RowVector3d& Rb,
                                                  const Eigen::Vector3d& Rab,
                                                  const Eigen::Vector3d& dipoleEvaluationCoordinate) {
  auto const& exponentGTFonA = gtoA.getExponent(GaussTypeFunctionIndexA);
  auto const& exponentGTFonB = gtoB.getExponent(GaussTypeFunctionIndexB);
  auto const exponentSum = exponentGTFonA + exponentGTFonB;
  initialize(Ra, Rb, dipoleEvaluationCoordinate, exponentGTFonA, exponentGTFonB, gtoA, gtoB);
  double exponentialCoefficientS_00 = -exponentGTFonB * exponentGTFonA / exponentSum;

  auto const normalizedCoefficientPrefactor = getNormalizationFactorAndCoefficient(
      GaussTypeFunctionIndexA, GaussTypeFunctionIndexB, gtoA, gtoB, exponentialCoefficientS_00, Rab);

  calculateSingleGTFPairIntegralOverShell(gtoA, gtoB, exponentSum);

  if (method_ == IntegralMethod::ClosedForm) {
    calculateAnalyticalDipoleIntegrals(GaussTypeFunctionIndexA, GaussTypeFunctionIndexB, gtoA, gtoB, Ra, Rb,
                                       dipoleEvaluationCoordinate);
  }
  createBlockForOneGTFOverShell(normalizedCoefficientPrefactor);
}

void GTODipoleMatrixBlock::calculateSingleGTFPairIntegralOverShell(const Utils::GtoExpansion& gtoA,
                                                                   const Utils::GtoExpansion& gtoB, double exponentSum) {
  // There are two relations for the Obara-Saika construction of S values:
  // read in T. Helgaker, P. Jorgensen, J. Olsen, Molecular electronic-structure theory, 2012
  // S_i,j = X_PA * S_i-1,j + 1/(2expSum) * ((i-1)*S_i-2,j + j*S_i-1,j-1)
  // S_i,j = X_PB * S_i,j-1 + 1/(2expSum) * (i*S_i-1,j-1 + (j-1)*S_i,j-2)
  // D_i,j = X_PA * D_i-1,j + 1/(2expSum) * ((i-1)*D_i-2,j + j*D_i-1,j-1 + S_i-1,j)
  // D_i,j = X_PB * D_i,j-1 + 1/(2expSum) * (i*D_i-1,j-1 + (j-1)*D_i,j-2 + S_i,j-1)

  const double p = 1 / (2 * exponentSum);
  for (int dimension = 0; dimension < 3; ++dimension) {
    // calculate directly overlap and dipole matrix elements, as they all will be used.
    for (int i = 0; i <= gtoA.angularMomentum(); ++i) {
      for (int j = 0; j <= gtoB.angularMomentum(); ++j) {
        if (i != 0 || j != 0) {
          if (j == 0) { // first formulae
            overlapIntegral_[dimension][i][j] = PminusA_[dimension] * overlapIntegral_[dimension][i - 1][j];
            dipoleIntegral_[dimension][i][j] =
                PminusA_[dimension] * dipoleIntegral_[dimension][i - 1][j] + p * overlapIntegral_[dimension][i - 1][j];
            if (i > 1) {
              overlapIntegral_[dimension][i][j] += (i - 1.0) * p * overlapIntegral_[dimension][i - 2][j];
              dipoleIntegral_[dimension][i][j] += (i - 1.0) * p * dipoleIntegral_[dimension][i - 2][j];
            }
          }
          else { // second formulae
            overlapIntegral_[dimension][i][j] = PminusB_[dimension] * overlapIntegral_[dimension][i][j - 1];
            dipoleIntegral_[dimension][i][j] =
                PminusB_[dimension] * dipoleIntegral_[dimension][i][j - 1] + p * overlapIntegral_[dimension][i][j - 1];
            if (i > 0) {
              overlapIntegral_[dimension][i][j] += i * p * overlapIntegral_[dimension][i - 1][j - 1];
              dipoleIntegral_[dimension][i][j] += i * p * dipoleIntegral_[dimension][i - 1][j - 1];
            }
            if (j > 1) {
              overlapIntegral_[dimension][i][j] += (j - 1.0) * p * overlapIntegral_[dimension][i][j - 2];
              dipoleIntegral_[dimension][i][j] += (j - 1.0) * p * dipoleIntegral_[dimension][i][j - 2];
            }
          }
        }
      } // end j loop
    }   // end i loop
  }     // end dimension loop
}

double GTODipoleMatrixBlock::getNormalizationFactorAndCoefficient(int GaussTypeFunctionIndexA, int GaussTypeFunctionIndexB,
                                                                  const Utils::GtoExpansion& gtoA,
                                                                  const Utils::GtoExpansion& gtoB,
                                                                  double expCoefficientS_00, const Eigen::Vector3d& Rab) {
  double x_ab = Rab(0);
  double y_ab = Rab(1);
  double z_ab = Rab(2);

  auto exponentialPrefactor = exp(expCoefficientS_00 * (x_ab * x_ab + y_ab * y_ab + z_ab * z_ab));
  auto GTFNormalizationCoefficients =
      gtoA.getNormalizedCoefficient(GaussTypeFunctionIndexA) * gtoB.getNormalizedCoefficient(GaussTypeFunctionIndexB);
  return exponentialPrefactor * GTFNormalizationCoefficients;
}

void GTODipoleMatrixBlock::createBlockForOneGTFOverShell(double normalizedCoefficientPrefactor) {
  for (int dipoleComponent = 0; dipoleComponent < 3; ++dipoleComponent) {
    for (int orbitalOnA = 0; orbitalOnA < orbitalShellSizeA_; ++orbitalOnA) {
      auto const AOIndexOnA = startAOonA_ + orbitalOnA;
      auto const AOMomentumOnA = AOMomenta_[AOIndexOnA];
      for (int orbitalOnB = 0; orbitalOnB < orbitalShellSizeB_; ++orbitalOnB) {
        auto const AOIndexOnB = startAOonB_ + orbitalOnB;
        auto const AOMomentumOnB = AOMomenta_[AOIndexOnB];
        auto X = overlapIntegral_[0][AOMomentumOnA.x][AOMomentumOnB.x];
        auto Y = overlapIntegral_[1][AOMomentumOnA.y][AOMomentumOnB.y];
        auto Z = overlapIntegral_[2][AOMomentumOnA.z][AOMomentumOnB.z];

        if (dipoleComponent == 0)
          X = dipoleIntegral_[0][AOMomentumOnA.x][AOMomentumOnB.x];
        else if (dipoleComponent == 1)
          Y = dipoleIntegral_[1][AOMomentumOnA.y][AOMomentumOnB.y];
        else
          Z = dipoleIntegral_[2][AOMomentumOnA.z][AOMomentumOnB.z];

        double const summationElement = normalizedCoefficientPrefactor * X * Y * Z;
        dipoleComponentsBlocks_[dipoleComponent](orbitalOnA, orbitalOnB) += summationElement;
      }
    }
  }
}

void GTODipoleMatrixBlock::dOrbitalsFromSixCartesianToFiveRealSolidHarmonics() {
  for (int dipoleComponent = 0; dipoleComponent < 3; ++dipoleComponent) {
    auto& dipoleBlock = dipoleComponentsBlocks_[dipoleComponent];

    if (orbitalShellSizeA_ == 6) {
      for (int j = 0; j < orbitalShellSizeB_; ++j) {
        dipoleBlock(3, j) = (dipoleBlock(3, j) - 0.5 * dipoleBlock(4, j) - 0.5 * dipoleBlock(5, j)) / sqrt(3);
        dipoleBlock(4, j) = 0.5 * (dipoleBlock(5, j) - dipoleBlock(4, j));
      }
    }
    if (orbitalShellSizeB_ == 6) {
      for (int i = 0; i < orbitalShellSizeA_; ++i) {
        dipoleBlock(i, 3) = (dipoleBlock(i, 3) - 0.5 * dipoleBlock(i, 4) - 0.5 * dipoleBlock(i, 5)) / sqrt(3);
        dipoleBlock(i, 4) = 0.5 * (dipoleBlock(i, 5) - dipoleBlock(i, 4));
      }
    }
  }
}

std::array<Eigen::MatrixXd, 3>
GTODipoleMatrixBlock::createSTOBlock(const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB,
                                     const Eigen::RowVector3d& Ra, const Eigen::RowVector3d& Rb,
                                     const Eigen::Vector3d& Rab, const Eigen::Vector3d& dipoleEvaluationCoordinate) {
  auto const numberGTFsOnA = gtoA.size();
  auto const numberGTFsOnB = gtoB.size();
  auto const blockRowSize = gtoA.nAOs();
  auto const blockColSize = gtoB.nAOs();

  for (auto& dCB : dipoleComponentsBlocks_)
    dCB.setConstant(6, 6, 0);

  for (int GTFonA = 0; GTFonA < numberGTFsOnA; ++GTFonA) {
    for (int GTFonB = 0; GTFonB < numberGTFsOnB; ++GTFonB) {
      calculateSingleGTFPair(GTFonA, GTFonB, gtoA, gtoB, Ra, Rb, Rab, dipoleEvaluationCoordinate);
    }
  }

  dOrbitalsFromSixCartesianToFiveRealSolidHarmonics();

  std::array<Eigen::MatrixXd, 3> dipoleBlocks = {};

  for (int dimension = 0; dimension < 3; ++dimension) {
    dipoleBlocks[dimension].resize(blockRowSize, blockColSize);
    dipoleBlocks[dimension] = dipoleComponentsBlocks_[dimension].topLeftCorner(blockRowSize, blockColSize);
  }

  return dipoleBlocks;
}

void GTODipoleMatrixBlock::calculateAnalyticalDipoleIntegrals(int GaussTypeFunctionIndexA, int GaussTypeFunctionIndexB,
                                                              const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB,
                                                              const Eigen::RowVector3d& Ra, const Eigen::RowVector3d& Rb,
                                                              const Eigen::Vector3d& evaluationCoordinates) {
  const double expA = gtoA.getExponent(GaussTypeFunctionIndexA);
  const double expB = gtoB.getExponent(GaussTypeFunctionIndexB);
  const double expSum = expA + expB;
  const double exponentialPrefactor = std::sqrt(M_PI / expSum) / expSum;
  twoDimensionalArray xDipoleIntegral, yDipoleIntegral, zDipoleIntegral;
  for (int i = 0; i <= gtoA.angularMomentum(); ++i) {
    for (int j = 0; j <= gtoB.angularMomentum(); ++j) {
      AnalyticalDipoleIntegralOverGTOsCalculator dipoleMatrixCalculator(i, j, expA, expB, Ra, Rb, evaluationCoordinates);
      auto const dipoleArray = dipoleMatrixCalculator.calculateAnalyticalDipoleElement();
      for (int dimension = 0; dimension < 3; ++dimension) {
        dipoleIntegral_[dimension][i][j] = dipoleArray[dimension] * exponentialPrefactor;
      }
    }
  }
}

void GTODipoleMatrixBlock::setIntegralMethod(IntegralMethod method) {
  method_ = method;
}

void GTODipoleMatrixBlock::initialize(const Eigen::RowVector3d& Ra, const Eigen::RowVector3d& Rb,
                                      const Eigen::RowVector3d& dipoleEvaluationCoordinate, double expA, double expB,
                                      const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB) {
  auto const exponentSum = expA + expB;
  auto const recursionInitialOverlap = std::sqrt(M_PI / exponentSum);
  for (int dimension = 0; dimension < 3; ++dimension) {
    overlapIntegral_[dimension][0][0] = recursionInitialOverlap;
    if (method_ == IntegralMethod::ObaraSaika) {
      PminusC_[dimension] =
          (Ra[dimension] * expA + Rb[dimension] * expB) / exponentSum - dipoleEvaluationCoordinate(dimension);
      dipoleIntegral_[dimension][0][0] = PminusC_[dimension] * overlapIntegral_[dimension][0][0];
    }
    PminusA_[dimension] = ((Ra[dimension] * expA + Rb[dimension] * expB) / exponentSum) - Ra[dimension];
    PminusB_[dimension] = ((Ra[dimension] * expA + Rb[dimension] * expB) / exponentSum) - Rb[dimension];
  }
  startAOonA_ = (gtoA.nAOs() == 1) ? 0 : (gtoA.nAOs() == 3) ? 1 : 4;
  startAOonB_ = (gtoB.nAOs() == 1) ? 0 : (gtoB.nAOs() == 3) ? 1 : 4;
  orbitalShellSizeA_ = (gtoA.nAOs() == 5) ? 6 : gtoA.nAOs();
  orbitalShellSizeB_ = (gtoB.nAOs() == 5) ? 6 : gtoB.nAOs();
}

} // namespace Sparrow
} // namespace Scine
