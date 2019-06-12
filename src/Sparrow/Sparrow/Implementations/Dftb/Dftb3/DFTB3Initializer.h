/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB3INITIALIZER_H
#define SPARROW_DFTB3INITIALIZER_H

#include "Sparrow/Implementations/Dftb/Utils/DFTBTypeInitializer.h"

namespace Scine {

namespace Utils {
class SinglePointMethod;
}

namespace Sparrow {
namespace dftb {
class DFTB3;

class DFTB3Initializer : public DFTBTypeInitializer {
 public:
  explicit DFTB3Initializer(Utils::SinglePointMethod* method);

 protected:
  void initialize(const Utils::PositionCollection& positions, const Utils::ElementTypeCollection& elementTypes) override;
  void initialize(const Utils::ElementTypeCollection& elementTypes) override;

 private:
  DFTB3& dftb3_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB3INITIALIZER_H