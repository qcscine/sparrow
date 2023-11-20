/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_CISMATRIXAOFOCKBUILDER_H
#define SPARROW_CISMATRIXAOFOCKBUILDER_H

#include "CISPseudoDensityBuilder.h"
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/Math/IterativeDiagonalizer/SigmaVectorEvaluator.h>

namespace Scine {
namespace Sparrow {
class CISData;

template<Utils::Reference restrictedness>
class CISMatrixAOFockBuilderBase {
 public:
  virtual Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>
  getAOFock(const Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>& pseudoDensity,
            std::map<int, std::vector<int>> atomPairList) const = 0;
  virtual ~CISMatrixAOFockBuilderBase() = default;
};

template<Utils::Reference restrictedness, Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet>
class CISMatrixAOFockBuilder : public CISMatrixAOFockBuilderBase<restrictedness> {
 public:
  using SAMType = Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>;
  CISMatrixAOFockBuilder(CISData cisData, const ExcitedStatesParam& excitedStatesParam);
  ~CISMatrixAOFockBuilder();
  SAMType getAOFock(const SAMType& pseudoDensity, std::map<int, std::vector<int>> atomPairList) const final;

 private:
  SAMType buildFock(const SAMType& pseudoDensity, std::map<int, std::vector<int>> atomPairList) const;
  void initialize();
  void calculateMatrices();
  void calculate(int atomI, int atomJ);
  void calculate(int atomI);
  Eigen::MatrixXd getCoulombIntegrals(int atomI, int atomJ);
  Eigen::MatrixXd getExchangeIntegrals(int atomI, int atomJ);
  CISData cisData_;
  int nAtoms_{cisData_.AOInfo.getNAtoms()};
  int nAOs_{cisData_.AOInfo.getNAtomicOrbitals()};
  std::shared_ptr<std::vector<std::map<int, std::shared_ptr<Eigen::MatrixXd>>>> coulombContainer_, exchangeContainer_;
  static constexpr const double sparseThreshold_ = 1e-8;
  double c1_;
  double c2_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_CISMATRIXAOFOCKBUILDER_H
