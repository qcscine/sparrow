/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SecondOrderFock.h"
#include "ZeroOrderMatricesCalculator.h"
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>

namespace Scine {
namespace Sparrow {

using std::exp;
using namespace Utils::AutomaticDifferentiation;

namespace dftb {

SecondOrderFock::SecondOrderFock(ZeroOrderMatricesCalculator& matricesCalculator,
                                 const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                                 const DFTBCommon::AtomicParameterContainer& atomicPar,
                                 const DFTBCommon::DiatomicParameterContainer& diatomicPar,
                                 const Utils::DensityMatrix& densityMatrix,
                                 const Eigen::MatrixXd& energyWeightedDensityMatrix, std::vector<double>& atomicCharges,
                                 const std::vector<double>& coreCharges, const Utils::AtomsOrbitalsIndexes& aoIndexes,
                                 const Eigen::MatrixXd& overlapMatrix, const bool& unrestrictedCalculationRunning)
  : ScfFock(matricesCalculator, elements, positions, atomicPar, diatomicPar, densityMatrix, energyWeightedDensityMatrix,
            atomicCharges, coreCharges, aoIndexes, overlapMatrix, unrestrictedCalculationRunning) {
}

void SecondOrderFock::initialize() {
  ScfFock::initialize();
  auto numberAtoms = elements_.size();

  G = Eigen::MatrixXd::Zero(numberAtoms, numberAtoms);
  dG.setDimension(numberAtoms, numberAtoms);
}

void SecondOrderFock::constructG(Utils::derivOrder order) {
  if (order == Utils::derivOrder::zero)
    constructG<Utils::derivOrder::zero>();
  else if (order == Utils::derivOrder::one)
    constructG<Utils::derivOrder::one>();
  else if (order == Utils::derivOrder::two)
    constructG<Utils::derivOrder::two>();
}

template<Utils::derivOrder O>
void SecondOrderFock::constructG() {
  // Constructs the G matrix, which contains the different
  // gamma_ab values for the different atom pairs, as well
  // as its derivative
  dG.setOrder(O);
  Value1DType<O> result;
#pragma omp parallel for private(result)
  for (int a = 0; a < getNumberAtoms(); ++a) {
    for (int b = a; b < getNumberAtoms(); ++b) {
      Eigen::Vector3d R = (positions_.row(b) - positions_.row(a));
      result = gamma<O>(a, b);
      auto v = get3Dfrom1D<O>(result, R);
#pragma omp critical(aFirst)
      { dG.get<O>()(a, b) = v; }
#pragma omp critical(bFirst)
      { dG.get<O>()(b, a) = getValueWithOppositeDerivative<O>(v); }
    }
  }
  G = dG.getMatrixXd();
}

template<Utils::derivOrder O>
Value1DType<O> SecondOrderFock::gamma(int a, int b) const {
  // Calculation of gamma according to elstner1998,
  // formulae are better explained in supplementary info of gaus2011
  auto R = variableWithUnitDerivative<O>((positions_.row(b) - positions_.row(a)).norm());
  auto R2 = R * R;
  double Ua = atomicPar_[Utils::ElementInfo::Z(elements_[a])]->getHubbardParameter();
  double Ub = atomicPar_[Utils::ElementInfo::Z(elements_[b])]->getHubbardParameter();

  if (getValue1DAsDouble<O>(R) == 0.0)
    return constant1D<O>(Ua);

  double ta = Ua * 3.2;
  double tb = Ub * 3.2;

  auto expa = exp(-ta * R);
  auto expb = exp(-tb * R);

  if (elements_[a] == elements_[b]) {
    double ta2 = ta * ta;
    auto expr = 1 / R + 0.6875 * ta + 0.1875 * R * ta2 + 0.02083333333333333 * R2 * ta * ta2; // From koehler2003, or
                                                                                              // supplementary info in
                                                                                              // gaus2011
    auto gamma = 1 / R - expa * expr;
    return gamma;
  }

  // Get precomputed parameters
  double term1a, term1b, term2a, term2b;
  diatomicPar_[Utils::ElementInfo::Z(elements_[a])][Utils::ElementInfo::Z(elements_[b])]->getGammaTerms(term1a, term1b,
                                                                                                        term2a, term2b);

  auto terma = -expa * (term1a - term2a / R);
  auto termb = -expb * (term1b - term2b / R);

  auto gamma = 1.0 / R + terma + termb;

  return gamma;
}

void SecondOrderFock::completeH() {
  correctionToFock.setZero();

  // The function completes H, meaning it calculates H = H0+H1
#pragma omp parallel for
  for (int a = 0; a < getNumberAtoms(); ++a) {
    int nAOsA = aoIndexes_.getNOrbitals(a);
    int indexA = aoIndexes_.getFirstOrbitalIndex(a);

    for (int b = a; b < getNumberAtoms(); b++) {
      int nAOsB = aoIndexes_.getNOrbitals(b);
      int indexB = aoIndexes_.getFirstOrbitalIndex(b);

      double sumOverAtoms = 0.0;
      for (int i = 0; i < getNumberAtoms(); i++)
        sumOverAtoms -= (G(a, i) + G(b, i)) * atomicCharges_[i];
      for (int mu = 0; mu < nAOsA; mu++) {
        for (int nu = 0; nu < nAOsB; nu++) {
#pragma omp atomic write
          HXoverS_(indexA + mu, indexB + nu) = 0.5 * sumOverAtoms;
#pragma omp atomic write
          correctionToFock(indexA + mu, indexB + nu) =
              overlapMatrix_(indexA + mu, indexB + nu) * HXoverS_(indexA + mu, indexB + nu);
          if (indexA != indexB) {
#pragma omp atomic read
            HXoverS_(indexB + nu, indexA + mu) = HXoverS_(indexA + mu, indexB + nu);
#pragma omp atomic read
            correctionToFock(indexB + nu, indexA + mu) = correctionToFock(indexA + mu, indexB + nu);
          }
        }
      }
    }
  }
}

void SecondOrderFock::addDerivatives(DerivativeContainerType<Utils::derivativeType::first>& derivatives) const {
  zeroOrderMatricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_ -
                                                               HXoverS_.cwiseProduct(densityMatrix_.restrictedMatrix()));
  addSecondOrderDerivatives<Utils::derivativeType::first>(derivatives);
  ScfFock::addDerivatives(derivatives);
}

void SecondOrderFock::addDerivatives(DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const {
  zeroOrderMatricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_ -
                                                               HXoverS_.cwiseProduct(densityMatrix_.restrictedMatrix()));
  addSecondOrderDerivatives<Utils::derivativeType::second_atomic>(derivatives);
  ScfFock::addDerivatives(derivatives);
}

void SecondOrderFock::addDerivatives(DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const {
  zeroOrderMatricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_ -
                                                               HXoverS_.cwiseProduct(densityMatrix_.restrictedMatrix()));
  addSecondOrderDerivatives<Utils::derivativeType::second_full>(derivatives);
  ScfFock::addDerivatives(derivatives);
}

template<Utils::derivativeType O>
void SecondOrderFock::addSecondOrderDerivatives(DerivativeContainerType<O>& derivatives) const {
  Value3DType<UnderlyingOrder<O>> der;
  DerivativeType<O> derivative;
  derivative.setZero();
#pragma omp parallel for firstprivate(derivative) private(der)
  for (int a = 0; a < getNumberAtoms(); ++a) {
    for (int b = a + 1; b < getNumberAtoms(); b++) {
      der = atomicCharges_[a] * atomicCharges_[b] * dG.get<UnderlyingOrder<O>>()(a, b);
      derivative = getDerivativeFromValueWithDerivatives<O>(der);
#pragma omp critical
      { addDerivativeToContainer<O>(derivatives, a, b, derivative); }
    }
  }

  if (unrestrictedCalculationRunning_)
    spinDFTB.addDerivatives<O>(derivatives, zeroOrderMatricesCalculator_.getOverlap(), densityMatrix_.alphaMatrix(),
                               densityMatrix_.betaMatrix());
}

double SecondOrderFock::calculateElectronicEnergy() const {
  auto numberAtoms = elements_.size();
  double elEnergy = (H0_.cwiseProduct(densityMatrix_.restrictedMatrix())).sum();

  for (int a = 0; a < numberAtoms; a++) {
    elEnergy += 0.5 * atomicCharges_[a] * atomicCharges_[a] * G(a, a);
#pragma omp simd reduction(+ : elEnergy)
    for (int b = a + 1; b < numberAtoms; b++)
      elEnergy += atomicCharges_[a] * atomicCharges_[b] * G(a, b);
  }

  if (unrestrictedCalculationRunning_) {
    elEnergy += spinDFTB.spinEnergyContribution();
  }

  return elEnergy;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
