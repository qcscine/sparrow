/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Dftb0Specifier.h"
#include "DFTB0.h"
#include "DFTB0Initializer.h"

namespace Scine {
namespace Sparrow {

namespace dftb {

std::string Dftb0Specifier::getName() const {
  return "DFTB0";
}

bool Dftb0Specifier::isLcaoMethod() const {
  return true;
}

bool Dftb0Specifier::isScfMethod() const {
  return false;
}

bool Dftb0Specifier::hasOrthonormalBasis() const {
  return false;
}

unsigned int Dftb0Specifier::maximalDerivativeOrder() const {
  return 2;
}

bool Dftb0Specifier::compatibleType(Utils::SinglePointMethod* method) const {
  return Utils::checkMethodType<DFTB0>(method);
}

std::unique_ptr<Utils::SinglePointMethod> Dftb0Specifier::createMethod() const {
  return std::make_unique<DFTB0>();
}

std::unique_ptr<Utils::MethodInitializer> Dftb0Specifier::createInitializer(Utils::SinglePointMethod* method) const {
  return std::make_unique<DFTB0Initializer>(method);
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
