/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB2INITIALIZER_H
#define SPARROW_DFTB2INITIALIZER_H

#include "Sparrow/Implementations/Dftb/Utils/DFTBTypeInitializer.h"

namespace Scine {
namespace Utils {
class SinglePointMethod;
}
namespace Sparrow {

namespace dftb {
class DFTB2;

class DFTB2Initializer : public DFTBTypeInitializer {
 public:
  explicit DFTB2Initializer(Utils::SinglePointMethod* method);

 protected:
  void initialize(const Utils::PositionCollection& positions, const Utils::ElementTypeCollection& elementTypes) override;
  void initialize(const Utils::ElementTypeCollection& elementTypes) override;

 private:
  DFTB2& dftb2_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB2INITIALIZER_H