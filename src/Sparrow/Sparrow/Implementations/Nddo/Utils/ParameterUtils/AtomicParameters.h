/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MNDOTYPEATOMICPARAMETERS_H
#define SPARROW_MNDOTYPEATOMICPARAMETERS_H

#include "PrincipalQuantumNumbers.h"
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ChargeSeparationParameter.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/KlopmanParameter.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/Geometry/ElementTypes.h>
#include <cmath>
#include <tuple>
#include <utility>

namespace Scine {
namespace Sparrow {
namespace nddo {

/*!
 * Class for the storage of atomic parameters in the
 * semiempirical methods. (only those needed at runtime).
 */

class AtomicParameters {
 public:
  explicit AtomicParameters(Utils::ElementType e = Utils::ElementType::none,
                            BasisFunctions basisFunctions = BasisFunctions::spd) {
    setElement(e, basisFunctions);
  }

  /*! Returns false if the element of the instance is not set. */
  bool isValid() const {
    return e_ != Utils::ElementType::none;
  }

  void setElement(Utils::ElementType e, BasisFunctions basisFunctions = BasisFunctions::spd) {
    e_ = e;
    setNAOs(PM6Elements::getNumberOfAOs(e, basisFunctions));
    setCoreCharge(PM6Elements::getCoreCharge(e));
  }

  Utils::ElementType element() const {
    return e_;
  }

  void setGTOs(Utils::AtomicGtos gtos) {
    gtos_ = std::move(gtos);
  }

  const Utils::AtomicGtos& GTOs() const {
    return gtos_;
  }

  void setNAOs(int n) {
    nAOs_ = n;
  }

  int nAOs() const {
    return nAOs_;
  }

  void setCoreCharge(double c) {
    coreCharge_ = c;
    coreCharge13_ = std::pow(c, 1. / 3.);
  }

  double coreCharge() const {
    return coreCharge_;
  }

  double cubicRootOfCoreCharge() const {
    return coreCharge13_;
  }

  void setBetaS(double v) {
    betaS_ = v;
  }

  void setBetaP(double v) {
    betaP_ = v;
  }

  void setBetaD(double v) {
    betaD_ = v;
  }

  double betaS() const {
    return betaS_;
  }

  double betaP() const {
    return betaP_;
  }

  double betaD() const {
    return betaD_;
  }

  void setUss(double v) {
    Uss_ = v;
  }

  void setUpp(double v) {
    Upp_ = v;
  }

  void setUdd(double v) {
    Udd_ = v;
  }

  double Uss() const {
    return Uss_;
  }

  double Upp() const {
    return Upp_;
  }

  double Udd() const {
    return Udd_;
  }

  double alpha() const {
    return alpha_;
  }

  void setAlpha(double v) {
    alpha_ = v;
  }

  void setPCore(double p) {
    pCore_ = p;
  }

  double pCore() const {
    return pCore_;
  }

  bool pCoreSpecified() const {
    return pCoreSpecified_;
  }

  void setPCoreSpecified(bool b) {
    pCoreSpecified_ = b;
  }

  void setKlopmanParameters(const multipole::KlopmanParameter& k) {
    klopman_ = k;
  }

  const multipole::KlopmanParameter& klopmanParameters() const {
    return klopman_;
  }

  void setChargeSeparations(const multipole::ChargeSeparationParameter& d) {
    chargeSep_ = d;
  }

  const multipole::ChargeSeparationParameter& chargeSeparations() const {
    return chargeSep_;
  }

  void addGaussianRepulsionParameters(double a, double b, double c) {
    gaussianRepulsionParameters_.emplace_back(a, b, c);
    hasGaussianRepulsionParameters_ = true;
  }

  bool hasGaussianRepulsionParameters() const {
    return hasGaussianRepulsionParameters_;
  }

  const std::vector<std::tuple<double, double, double>>& getGaussianRepulsionParameters() const {
    return gaussianRepulsionParameters_;
  }

  void clearGaussianRepulsionParameters() {
    gaussianRepulsionParameters_.clear();
  }

 private:
  Utils::ElementType e_{Utils::ElementType::none};
  int nAOs_{0};
  double coreCharge_{0.};
  double coreCharge13_{}; // cubic root of coreCharge_
  Utils::AtomicGtos gtos_{};
  double betaS_{0.}, betaP_{0.}, betaD_{0.};
  double Uss_{0.}, Upp_{0.}, Udd_{0.};
  double alpha_{0.};
  double pCore_{0.};
  bool hasGaussianRepulsionParameters_{false};
  std::vector<std::tuple<double, double, double>> gaussianRepulsionParameters_{}; // NB: make sure they have units
                                                                                  // compatible with formula so that R
                                                                                  // does not need to be converted to
                                                                                  //
  // std::tuple<double, double, double> gaussianRepulsionParameters_; // NB: make sure they have units compatible with
  // formula so that R does not need to be converted to
  // another unit!
  multipole::KlopmanParameter klopman_{};
  multipole::ChargeSeparationParameter chargeSep_{};
  bool pCoreSpecified_{false};
};

} // namespace nddo
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_MNDOTYPEATOMICPARAMETERS_H
