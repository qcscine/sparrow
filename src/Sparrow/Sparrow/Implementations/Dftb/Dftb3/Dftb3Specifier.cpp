/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Dftb3Specifier.h"
#include "DFTB3.h"
#include "DFTB3Initializer.h"

namespace Scine {
namespace Sparrow {

namespace dftb {

std::string Dftb3Specifier::getName() const {
  return "DFTB3";
}

bool Dftb3Specifier::isLcaoMethod() const {
  return true;
}

bool Dftb3Specifier::isScfMethod() const {
  return true;
}

bool Dftb3Specifier::hasOrthonormalBasis() const {
  return false;
}

unsigned int Dftb3Specifier::maximalDerivativeOrder() const {
  return 2;
}

bool Dftb3Specifier::compatibleType(Utils::SinglePointMethod* method) const {
  return Utils::checkMethodType<DFTB3>(method);
}

std::unique_ptr<Utils::SinglePointMethod> Dftb3Specifier::createMethod() const {
  return std::make_unique<DFTB3>();
}

std::unique_ptr<Utils::MethodInitializer> Dftb3Specifier::createInitializer(Utils::SinglePointMethod* method) const {
  return std::make_unique<DFTB3Initializer>(method);
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
