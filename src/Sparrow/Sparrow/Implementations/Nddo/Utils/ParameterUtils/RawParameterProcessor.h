/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_RAWPARAMETERPROCESSOR_H
#define SPARROW_RAWPARAMETERPROCESSOR_H

#include "PrincipalQuantumNumbers.h"
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/SlaterCondonParameters.h>
#include <Utils/Geometry/ElementTypes.h>
#include <memory>

namespace Scine {
namespace Sparrow {
namespace nddo {

class RawAtomicParameters;
class RawDiatomicParameters;
class RawParametersContainer;
class AtomicParameters;
class PM6DiatomicParameters;
class OneCenterTwoElectronIntegrals;
/*!
 * This class implements functions for the conversion between
 * raw parameters published for PM6 and parameters useful at
 * runtime.
 */
class RawParameterProcessor {
 public:
  explicit RawParameterProcessor(const RawParametersContainer& rawParameters,
                                 BasisFunctions basisFunctions = BasisFunctions::spd);
  std::unique_ptr<PM6DiatomicParameters> runtimeDiatomicParameters(Utils::ElementType e1, Utils::ElementType e2);
  std::pair<std::unique_ptr<AtomicParameters>, std::unique_ptr<OneCenterTwoElectronIntegrals>>
  processAtomicParameters(Utils::ElementType e);

 private:
  std::unique_ptr<OneCenterTwoElectronIntegrals> get1c2eIntegrals(Utils::ElementType e, const RawAtomicParameters& p) const;
  void computeSlaterCondonParameters(AtomicParameters& runtimeAtomicPar, const RawAtomicParameters& p);
  void setKlopman(AtomicParameters& par, const RawAtomicParameters& p) const;
  void setChargeSeparations(Utils::ElementType e, AtomicParameters& par, const RawAtomicParameters& p) const;
  void setGtoExpansion(Utils::ElementType e, AtomicParameters& par, const RawAtomicParameters& p) const;
  void setDiatomicExponent(PM6DiatomicParameters& par, Utils::ElementType e1, Utils::ElementType e2,
                           const RawDiatomicParameters& p);
  const RawParametersContainer& rawParameters_;
  SlaterCondonParameters scParameters_;
  BasisFunctions basisFunctions_;
};
} // namespace nddo
} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_RAWPARAMETERPROCESSOR_H
