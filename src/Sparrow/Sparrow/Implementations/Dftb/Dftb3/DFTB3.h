/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB3_H
#define SPARROW_DFTB3_H

#include "Sparrow/Implementations/Dftb/Utils/DFTBCommon.h"
#include "Sparrow/Implementations/Dftb/Utils/SDFTB.h"
#include "Utils/DataStructures/MatrixWithDerivatives.h"
#include <Utils/Scf/MethodInterfaces/ScfMethod.h>
#include <Eigen/Core>
#include <string>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace dftb {
class ZeroOrderMatricesCalculator;

class DFTB3 : public Utils::ScfMethod {
 public:
  DFTB3();
  ~DFTB3() override;
  void initializeFromParameterPath(const std::string& path);
  std::shared_ptr<DFTBCommon> getInitializer() const;
  Eigen::MatrixXd calculateGammaMatrix() const;
  std::shared_ptr<Eigen::VectorXd> calculateSpinConstantVector() const;

 private:
  DFTBCommon::AtomicParameterContainer atomParameters;   // parameters for atoms
  DFTBCommon::DiatomicParameterContainer pairParameters; // List of pointers to parameters
  std::shared_ptr<DFTBCommon> dftbBase;
  std::unique_ptr<dftb::ZeroOrderMatricesCalculator> matricesCalculator_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB3_H
