/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_RAWPARAMETERPROCESSOR_H
#define SPARROW_RAWPARAMETERPROCESSOR_H

#include "PrincipalQuantumNumbers.h"
#include <Sparrow/Implementations/Nddo/Parameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/SlaterCondonParameters.h>
#include <Utils/Geometry/ElementTypes.h>
#include <memory>

namespace Scine {
namespace Sparrow {
namespace nddo {

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
  explicit RawParameterProcessor(const Parameters& rawParameters, BasisFunctions basisFunctions = BasisFunctions::spd);
  std::unique_ptr<PM6DiatomicParameters> runtimeDiatomicParameters(Utils::ElementType e1, Utils::ElementType e2);
  std::pair<std::unique_ptr<AtomicParameters>, std::unique_ptr<OneCenterTwoElectronIntegrals>>
  processAtomicParameters(Utils::ElementType e);

 private:
  std::unique_ptr<OneCenterTwoElectronIntegrals> get1c2eIntegrals(Utils::ElementType e, const Parameters::Atomic& p) const;
  void computeSlaterCondonParameters(AtomicParameters& runtimeAtomicPar, const Parameters::Atomic& p);
  void setKlopman(AtomicParameters& par, const Parameters::Atomic& p) const;
  void setChargeSeparations(Utils::ElementType e, AtomicParameters& par, const Parameters::Atomic& p) const;
  void setGtoExpansion(Utils::ElementType e, AtomicParameters& par, const Parameters::Atomic& p) const;
  static void setDiatomicExponent(PM6DiatomicParameters& par, Utils::ElementType e1, Utils::ElementType e2,
                                  const Parameters::Diatomic& p);

  const Parameters& rawParameters_;
  SlaterCondonParameters scParameters_;
  BasisFunctions basisFunctions_;
};
} // namespace nddo
} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_RAWPARAMETERPROCESSOR_H
