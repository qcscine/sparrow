/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AtomPairDipole.h"
#include "GTODipoleMatrixBlock.h"
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMatrixCalculator.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/DipoleMatrix.h>

namespace Scine {
namespace Sparrow {

void AtomPairDipole::fillAtomPairDipoleBlock(Utils::DipoleMatrix& dipoleMatrix, int startOfAtomA, int startOfAtomB,
                                             const IntegralMethod& method, const Utils::AtomicGtos& gtosA,
                                             const Utils::AtomicGtos& gtosB, const Eigen::RowVector3d& Ra,
                                             const Eigen::RowVector3d& Rb, const Eigen::RowVector3d& Rab,
                                             const Eigen::RowVector3d& dipoleEvaluationCoordinate) {
  GTODipoleMatrixBlock block;
  block.setIntegralMethod(method);

  using Utils::DerivativeOrder;

  if (gtosA.s) {
    if (gtosB.s) {
      const auto ssBlock = block.createSTOBlock(gtosA.s.value(), gtosB.s.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB, 1, 1) = ssBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB, 1, 1) = ssBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB, 1, 1) = ssBlock[2];
    }
    if (gtosB.p) {
      const auto spBlock = block.createSTOBlock(gtosA.s.value(), gtosB.p.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB + 1, 1, 3) = spBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB + 1, 1, 3) = spBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB + 1, 1, 3) = spBlock[2];
    }
    if (gtosB.d) {
      const auto sdBlock = block.createSTOBlock(gtosA.s.value(), gtosB.d.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB + 4, 1, 5) = sdBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB + 4, 1, 5) = sdBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA, startOfAtomB + 4, 1, 5) = sdBlock[2];
    }
  }
  if (gtosA.p) {
    if (gtosB.s) {
      const auto psBlock = block.createSTOBlock(gtosA.p.value(), gtosB.s.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB, 3, 1) = psBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB, 3, 1) = psBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB, 3, 1) = psBlock[2];
    }
    if (gtosB.p) {
      const auto ppBlock = block.createSTOBlock(gtosA.p.value(), gtosB.p.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB + 1, 3, 3) = ppBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB + 1, 3, 3) = ppBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB + 1, 3, 3) = ppBlock[2];
    }
    if (gtosB.d) {
      const auto pdBlock = block.createSTOBlock(gtosA.p.value(), gtosB.d.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB + 4, 3, 5) = pdBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB + 4, 3, 5) = pdBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA + 1, startOfAtomB + 4, 3, 5) = pdBlock[2];
    }
  }
  if (gtosA.d) {
    if (gtosB.s) {
      const auto dsBlock = block.createSTOBlock(gtosA.d.value(), gtosB.s.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB, 5, 1) = dsBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB, 5, 1) = dsBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB, 5, 1) = dsBlock[2];
    }
    if (gtosB.p) {
      const auto dpBlock = block.createSTOBlock(gtosA.d.value(), gtosB.p.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB + 1, 5, 3) = dpBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB + 1, 5, 3) = dpBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB + 1, 5, 3) = dpBlock[2];
    }
    if (gtosB.d) {
      const auto ddBlock = block.createSTOBlock(gtosA.d.value(), gtosB.d.value(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB + 4, 5, 5) = ddBlock[0];
      dipoleMatrix.y().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB + 4, 5, 5) = ddBlock[1];
      dipoleMatrix.z().get<DerivativeOrder::Zero>().block(startOfAtomA + 4, startOfAtomB + 4, 5, 5) = ddBlock[2];
    }
  }
}
} // namespace Sparrow
} // namespace Scine
