/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ThirdOrderFock.h"
#include "Utils/DataStructures/DensityMatrix.h"
#include "Utils/Scf/LcaoUtils/LcaoUtils.h"
#include "ZeroOrderMatricesCalculator.h"
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>

namespace Scine {
namespace Sparrow {

using std::exp;
using namespace Utils::AutomaticDifferentiation;

namespace dftb {

ThirdOrderFock::ThirdOrderFock(ZeroOrderMatricesCalculator& matricesCalculator,
                               const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                               const DFTBCommon::AtomicParameterContainer& atomicPar,
                               const DFTBCommon::DiatomicParameterContainer& diatomicPar,
                               const Utils::DensityMatrix& densityMatrix,
                               const Eigen::MatrixXd& energyWeightedDensityMatrix, std::vector<double>& atomicCharges,
                               const std::vector<double>& coreCharges, const Utils::AtomsOrbitalsIndexes& aoIndexes,
                               const Eigen::MatrixXd& overlapMatrix, const bool& unrestrictedCalculationRunning)
  : ScfFock(matricesCalculator, elements, positions, atomicPar, diatomicPar, densityMatrix, energyWeightedDensityMatrix,
            atomicCharges, coreCharges, aoIndexes, overlapMatrix, unrestrictedCalculationRunning),
    zeta(4.0) {
}

void ThirdOrderFock::initialize() {
  ScfFock::initialize();
  auto numberAtoms = elements_.size();

  g = Eigen::MatrixXd::Zero(numberAtoms, numberAtoms);
  G = Eigen::MatrixXd::Zero(numberAtoms, numberAtoms);
  dg.setDimension(numberAtoms, numberAtoms);
  dG.setDimension(numberAtoms, numberAtoms);
}

void ThirdOrderFock::constructG(Utils::derivOrder order) {
  if (order == Utils::derivOrder::zero)
    constructG<Utils::derivOrder::zero>();
  else if (order == Utils::derivOrder::one)
    constructG<Utils::derivOrder::one>();
  else if (order == Utils::derivOrder::two)
    constructG<Utils::derivOrder::two>();
}

template<Utils::derivOrder O>
void ThirdOrderFock::constructG() {
  dg.setOrder(O);
  dG.setOrder(O);
  Value1DType<O> gamma, Gab, Gba;
#pragma omp parallel for private(gamma, Gab, Gba)
  for (int a = 0; a < getNumberAtoms(); ++a) {
    for (int b = a; b < getNumberAtoms(); b++) {
      Eigen::Vector3d R = (positions_.row(b) - positions_.row(a));
      gammah<O>(a, b, gamma, Gab, Gba);
      auto term = get3Dfrom1D<O>(gamma, R);
#pragma omp critical(dg)
      {
        dg.get<O>()(a, b) = term;
        dg.get<O>()(b, a) = getValueWithOppositeDerivative<O>(term);
      }
#pragma omp critical(dG)
      {
        dG.get<O>()(a, b) = get3Dfrom1D<O>(Gab, R);
        dG.get<O>()(b, a) = get3Dfrom1D<O>(Gba, R); // Not so sure...
      }
    }
  }
  g = dg.getMatrixXd();
  G = dG.getMatrixXd();
}

template<Utils::derivOrder O>
void ThirdOrderFock::gammah(int a, int b, Value1DType<O>& gamma, Value1DType<O>& Gab, Value1DType<O>& Gba) const {
  // Calculate gamma^h according to gaus2011
  // Notation is more or less the same as in its supplementary information
  auto R = variableWithUnitDerivative<O>((positions_.row(b) - positions_.row(a)).norm());
  auto R2 = R * R;
  double Ua = atomicPar_[Utils::ElementInfo::Z(elements_[a])]->getHubbardParameter();
  double Ub = atomicPar_[Utils::ElementInfo::Z(elements_[b])]->getHubbardParameter();

  if (aoIndexes_.getFirstOrbitalIndex(a) == aoIndexes_.getFirstOrbitalIndex(b)) {
    gamma = constant1D<O>(Ua);
    double v = 0.5 * atomicPar_[Utils::ElementInfo::Z(elements_[a])]->getHubbardDerivative();
    Gab = Gba = constant1D<O>(v);
    return;
  }

  // Calculate h and its derivatives
  Value1DType<O> h, dhdU;
  if (elements_[a] == Utils::ElementType::H || elements_[b] == Utils::ElementType::H)
    hFactor<O>(Ua, Ub, R, h, dhdU);
  else {
    h = constant1D<O>(1.0);
    dhdU = constant1D<O>(0.0);
  }

  // Get values that are employed several times
  double ta = Ua * 3.2;
  double tb = Ub * 3.2;
  auto expa = exp(-ta * R);
  auto expb = exp(-tb * R);

  // Calculate S = {Sg or Sf} and its derivatives
  Value1DType<O> S, dSda, dSdb;
  if (elements_[a] == elements_[b]) {
    double ta2 = ta * ta;
    // double tb2 = tb*tb; //never used
    auto g = 1.0 / 48.0 * (48.0 / R + 33.0 * ta + 9.0 * R * ta2 + R2 * ta * ta2); // From koehler2003, or supplementary
                                                                                  // info in gaus2011
    auto dgda = 1.0 / 48.0 * (33.0 + 18.0 * ta * R + 3.0 * ta2 * R2);

    S = expa * g;
    dSda = dSdb = expa * (dgda - R * g);
  }
  else {
    // Get precomputed values
    double term1a, term1b, term2a, term2b;
    diatomicPar_[Utils::ElementInfo::Z(elements_[a])][Utils::ElementInfo::Z(elements_[b])]->getGammaTerms(term1a, term1b,
                                                                                                          term2a, term2b);
    double dtermab1a, dtermab1b, dtermab2a, dtermab2b, dtermba1a, dtermba1b, dtermba2a, dtermba2b, dtadr, dtbdr;
    diatomicPar_[Utils::ElementInfo::Z(elements_[a])][Utils::ElementInfo::Z(elements_[b])]->getGammaDerTerms(
        dtermab1a, dtermab1b, dtermab2a, dtermab2b, dtermba1a, dtermba1b, dtermba2a, dtermba2b, dtadr, dtbdr);

    // Calculate derived values
    auto fab = (term1a - term2a / R);
    auto fba = (term1b - term2b / R);
    auto dfabda = -dtermab1a - dtermab2a / R;
    auto dfbadb = -dtermba1a - dtermba2a / R;
    auto dfabdb = dtermab1b + dtermab2b / R;
    auto dfbada = dtermba1b + dtermba2b / R;

    S = expa * fab + expb * fba;
    dSda = (expa * (dfabda - R * fab) + expb * dfbada);
    dSdb = (expb * (dfbadb - R * fba) + expa * dfabdb);
  }

  // Gamma derivatives:
  auto dgdUa = -(3.2 * dSda * h + S * dhdU);
  auto dgdUb = -(3.2 * dSdb * h + S * dhdU);

  // Set values to return
  gamma = 1.0 / R - S * h;
  Gab = dgdUa * atomicPar_[Utils::ElementInfo::Z(elements_[a])]->getHubbardDerivative();
  Gba = dgdUb * atomicPar_[Utils::ElementInfo::Z(elements_[b])]->getHubbardDerivative();
}

template<Utils::derivOrder O>
void ThirdOrderFock::hFactor(double Ua, double Ub, const Value1DType<O>& R, Value1DType<O>& h, Value1DType<O>& dhdU) const {
  // According to gaus2011
  auto R2 = R * R;
  auto powerM1 = std::pow((Ua + Ub) / 2.0, zeta - 1.0); // Calculate it only once
  auto power = powerM1 * (Ua + Ub) / 2.0;

  h = exp(-power * R2);
  dhdU = -zeta * R2 / 2.0 * powerM1 * h;
}

void ThirdOrderFock::completeH() {
  correctionToFock.setZero();

#pragma omp parallel for
  for (int a = 0; a < getNumberAtoms(); ++a) {
    int nAOsA = aoIndexes_.getNOrbitals(a);
    int indexA = aoIndexes_.getFirstOrbitalIndex(a);

    for (int b = a; b < getNumberAtoms(); b++) {
      int nAOsB = aoIndexes_.getNOrbitals(b);
      int indexB = aoIndexes_.getFirstOrbitalIndex(b);

      double H1sumOverAtoms = 0.0;
      double H2sumOverAtoms = 0.0;
      for (int i = 0; i < getNumberAtoms(); i++) {
        H1sumOverAtoms -= (g(a, i) + g(b, i)) * atomicCharges_[i];
        H2sumOverAtoms += atomicCharges_[i] * (2 * (atomicCharges_[a] * G(a, i) + atomicCharges_[b] * G(b, i)) +
                                               atomicCharges_[i] * (G(i, a) + G(i, b)));
      }
      for (int mu = 0; mu < nAOsA; mu++) {
        for (int nu = 0; nu < nAOsB; nu++) {
          HXoverS_(indexA + mu, indexB + nu) = 0.5 * H1sumOverAtoms + H2sumOverAtoms / 6.0;
          correctionToFock(indexA + mu, indexB + nu) =
              overlapMatrix_(indexA + mu, indexB + nu) * HXoverS_(indexA + mu, indexB + nu);
          if (indexA != indexB) {
            HXoverS_(indexB + nu, indexA + mu) = HXoverS_(indexA + mu, indexB + nu);
            correctionToFock(indexB + nu, indexA + mu) = correctionToFock(indexA + mu, indexB + nu);
          }
        }
      }
    }
  }
}

void ThirdOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const {
  zeroOrderMatricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_ -
                                                               HXoverS_.cwiseProduct(densityMatrix_.restrictedMatrix()));
  addThirdOrderDerivatives<Utils::derivativeType::first>(derivatives);
  ScfFock::addDerivatives(derivatives);
}

void ThirdOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const {
  zeroOrderMatricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_ -
                                                               HXoverS_.cwiseProduct(densityMatrix_.restrictedMatrix()));
  addThirdOrderDerivatives<Utils::derivativeType::second_atomic>(derivatives);
  ScfFock::addDerivatives(derivatives);
}

void ThirdOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const {
  zeroOrderMatricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_ -
                                                               HXoverS_.cwiseProduct(densityMatrix_.restrictedMatrix()));
  addThirdOrderDerivatives<Utils::derivativeType::second_full>(derivatives);
  ScfFock::addDerivatives(derivatives);
}

template<Utils::derivativeType O>
void ThirdOrderFock::addThirdOrderDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const {
  Value3DType<UnderlyingOrder<O>> der;
  DerivativeType<O> derivative;
  derivative.setZero();
#pragma omp parallel for firstprivate(derivative) private(der)
  for (int a = 0; a < getNumberAtoms(); ++a) {
    for (int b = a + 1; b < getNumberAtoms(); b++) {
      der = atomicCharges_[a] * atomicCharges_[b] * dg.get<UnderlyingOrder<O>>()(b, a); // nb not dg[a][b] because we
                                                                                        // want dgab/dRa
      der += -1. * atomicCharges_[a] / 3.0 * atomicCharges_[b] *
             (atomicCharges_[b] * dG.get<UnderlyingOrder<O>>()(b, a) - atomicCharges_[a] * dG.get<UnderlyingOrder<O>>()(a, b));
      derivative = getDerivativeFromValueWithDerivatives<O>(der);
#pragma omp critical
      { addDerivativeToContainer<O>(derivatives, b, a, derivative); }
    }
  }

  if (unrestrictedCalculationRunning_)
    spinDFTB.addDerivatives<O>(derivatives, zeroOrderMatricesCalculator_.getOverlap(), densityMatrix_.alphaMatrix(),
                               densityMatrix_.betaMatrix());
}

double ThirdOrderFock::calculateElectronicEnergy() const {
  double elEnergy = (H0_.cwiseProduct(densityMatrix_.restrictedMatrix())).sum();

#pragma omp simd collapse(2) reduction(+ : elEnergy)
  for (int a = 0; a < getNumberAtoms(); a++) {
    for (int b = 0; b < getNumberAtoms(); b++) {
      elEnergy += 0.5 * atomicCharges_[a] * atomicCharges_[b] * g(a, b);
      elEnergy -= atomicCharges_[a] * atomicCharges_[a] * atomicCharges_[b] * G(a, b) / 3.0;
    }
  }

  if (unrestrictedCalculationRunning_) {
    elEnergy += spinDFTB.spinEnergyContribution();
  }

  return elEnergy;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
