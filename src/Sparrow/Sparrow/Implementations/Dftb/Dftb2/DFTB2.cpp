/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTB2.h"
#include "Sparrow/Implementations/Dftb/Utils/DFTBCommon.h"
#include "Sparrow/Implementations/Dftb/Utils/DensityGuess.h"
#include "Sparrow/Implementations/Dftb/Utils/Overlap.h"
#include "Sparrow/Implementations/Dftb/Utils/Repulsion.h"
#include "Sparrow/Implementations/Dftb/Utils/SecondOrderFock.h"
#include "Sparrow/Implementations/Dftb/Utils/ZeroOrderMatricesCalculator.h"
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace dftb {

DFTB2::DFTB2() : ScfMethod(true, Utils::derivOrder::two), atomParameters(110) {
  dftbBase = std::make_shared<DFTBCommon>(elementTypes_, nElectrons_, molecularCharge_, atomParameters, pairParameters);
  matricesCalculator_ = std::make_unique<dftb::ZeroOrderMatricesCalculator>(elementTypes_, positions_, aoIndexes_,
                                                                            atomParameters, pairParameters, densityMatrix_);
  overlapCalculator_ = std::make_unique<dftb::Overlap>(*matricesCalculator_);
  electronicPart_ = std::make_unique<dftb::SecondOrderFock>(*matricesCalculator_, elementTypes_, positions_,
                                                            atomParameters, pairParameters, densityMatrix_,
                                                            energyWeightedDensityMatrix_, atomicCharges_, coreCharges_,
                                                            aoIndexes_, overlapMatrix_, unrestrictedCalculationRunning_);
  rep_ = std::make_unique<dftb::Repulsion>(elementTypes_, positions_, dftbBase->getPairParameters());
  densityMatrixGuess_ = std::make_unique<dftb::DensityGuess>(aoIndexes_, coreCharges_, nElectrons_);

  initializer_ = dftbBase;
}

DFTB2::~DFTB2() = default;

void DFTB2::initializeFromParameterPath(const std::string& path) {
  dftbBase->setMethodDetails(path, 2);
  ScfMethod::initialize();
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
