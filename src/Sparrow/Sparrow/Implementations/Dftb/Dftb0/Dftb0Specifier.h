/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB0SPECIFIER_H
#define SPARROW_DFTB0SPECIFIER_H

#include <Utils/MethodEssentials/MethodFactories/MethodSpecifier.h>
#include <Utils/MethodEssentials/Methods/SinglePointMethod.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

class Dftb0Specifier : public Utils::MethodSpecifier {
 public:
  std::string getName() const override;
  bool isLcaoMethod() const override;
  bool isScfMethod() const override;
  bool hasOrthonormalBasis() const override;
  unsigned int maximalDerivativeOrder() const override;
  bool compatibleType(Utils::SinglePointMethod* method) const override;
  std::unique_ptr<Utils::SinglePointMethod> createMethod() const override;
  std::unique_ptr<Utils::MethodInitializer> createInitializer(Utils::SinglePointMethod* method) const override;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB0SPECIFIER_H
