/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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

  using Utils::derivOrder;

  if (gtosA.hasS()) {
    if (gtosB.hasS()) {
      const auto ssBlock = block.createSTOBlock(gtosA.s(), gtosB.s(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB, 1, 1) = ssBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB, 1, 1) = ssBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB, 1, 1) = ssBlock[2];
    }
    if (gtosB.hasP()) {
      const auto spBlock = block.createSTOBlock(gtosA.s(), gtosB.p(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB + 1, 1, 3) = spBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB + 1, 1, 3) = spBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB + 1, 1, 3) = spBlock[2];
    }
    if (gtosB.hasD()) {
      const auto sdBlock = block.createSTOBlock(gtosA.s(), gtosB.d(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB + 4, 1, 5) = sdBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB + 4, 1, 5) = sdBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA, startOfAtomB + 4, 1, 5) = sdBlock[2];
    }
  }
  if (gtosA.hasP()) {
    if (gtosB.hasS()) {
      const auto psBlock = block.createSTOBlock(gtosA.p(), gtosB.s(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB, 3, 1) = psBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB, 3, 1) = psBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB, 3, 1) = psBlock[2];
    }
    if (gtosB.hasP()) {
      const auto ppBlock = block.createSTOBlock(gtosA.p(), gtosB.p(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB + 1, 3, 3) = ppBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB + 1, 3, 3) = ppBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB + 1, 3, 3) = ppBlock[2];
    }
    if (gtosB.hasD()) {
      const auto pdBlock = block.createSTOBlock(gtosA.p(), gtosB.d(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB + 4, 3, 5) = pdBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB + 4, 3, 5) = pdBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA + 1, startOfAtomB + 4, 3, 5) = pdBlock[2];
    }
  }
  if (gtosA.hasD()) {
    if (gtosB.hasS()) {
      const auto dsBlock = block.createSTOBlock(gtosA.d(), gtosB.s(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB, 5, 1) = dsBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB, 5, 1) = dsBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB, 5, 1) = dsBlock[2];
    }
    if (gtosB.hasP()) {
      const auto dpBlock = block.createSTOBlock(gtosA.d(), gtosB.p(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB + 1, 5, 3) = dpBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB + 1, 5, 3) = dpBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB + 1, 5, 3) = dpBlock[2];
    }
    if (gtosB.hasD()) {
      const auto ddBlock = block.createSTOBlock(gtosA.d(), gtosB.d(), Ra, Rb, Rab, dipoleEvaluationCoordinate);
      dipoleMatrix.x().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB + 4, 5, 5) = ddBlock[0];
      dipoleMatrix.y().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB + 4, 5, 5) = ddBlock[1];
      dipoleMatrix.z().get<derivOrder::zero>().block(startOfAtomA + 4, startOfAtomB + 4, 5, 5) = ddBlock[2];
    }
  }
}
} // namespace Sparrow
} // namespace Scine
