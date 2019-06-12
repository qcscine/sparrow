/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Dftb2Specifier.h"
#include "DFTB2.h"
#include "DFTB2Initializer.h"

namespace Scine {
namespace Sparrow {

namespace dftb {

std::string Dftb2Specifier::getName() const {
  return "DFTB2";
}

bool Dftb2Specifier::isLcaoMethod() const {
  return true;
}

bool Dftb2Specifier::isScfMethod() const {
  return true;
}

bool Dftb2Specifier::hasOrthonormalBasis() const {
  return false;
}

unsigned int Dftb2Specifier::maximalDerivativeOrder() const {
  return 2;
}

bool Dftb2Specifier::compatibleType(Utils::SinglePointMethod* method) const {
  return Utils::checkMethodType<DFTB2>(method);
}

std::unique_ptr<Utils::SinglePointMethod> Dftb2Specifier::createMethod() const {
  return std::make_unique<DFTB2>();
}

std::unique_ptr<Utils::MethodInitializer> Dftb2Specifier::createInitializer(Utils::SinglePointMethod* method) const {
  return std::make_unique<DFTB2Initializer>(method);
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
