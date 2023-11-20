/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal dependencies */
#include "TDDFTBCalculator.h"
#include "BasisPruner.h"
#include "TDDFTBData.h"
#include "TDDFTBEigenvalueSolver.h"
#include "TDDFTBSettings.h"
#include <Sparrow/Implementations/Dftb/DFTBMethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Dftb0/Wrapper/DFTB0MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/TransitionChargesCalculator.h>
#include <Sparrow/Implementations/Exceptions.h>
#include <Sparrow/Implementations/TimeDependent/DiagonalPreconditionerEvaluator.h>
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
/* External dependencies */
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/Constants.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/TimeDependent/TransitionDipoleCalculator.h>

namespace Scine {
namespace Sparrow {

InvalidCalculatorTypeForTDDFTB::InvalidCalculatorTypeForTDDFTB(std::shared_ptr<Scine::Core::Calculator> method) {
  error = "Calculator " + method->name() + " is not an dftb method and cannot calculate TDDFTB.";
}

void TDDFTBCalculator::setReferenceCalculator(std::shared_ptr<Scine::Core::Calculator> method) {
  dftbMethod_ = std::dynamic_pointer_cast<DFTBMethodWrapper>(method);
  if (!dftbMethod_) {
    throw InvalidCalculatorTypeForTDDFTB(method);
  }
  tddftbData_ = std::make_unique<TDDFTBData>(dftbMethod_->getTDDFTBData());
}

TDDFTBCalculator::TDDFTBCalculator() {
  settings_ = std::make_unique<TDDFTBSettings>();
}

Utils::Results& TDDFTBCalculator::results() {
  return results_;
}

const Utils::Results& TDDFTBCalculator::results() const {
  return results_;
}

void TDDFTBCalculator::referenceCalculation() {
  if (!dftbMethod_)
    throw MissingReferenceCalculatorException();
  if (dftbMethod_->settings().valueExists(Utils::SettingsNames::selfConsistenceCriterion) &&
      dftbMethod_->settings().getDouble(Utils::SettingsNames::selfConsistenceCriterion) > 1e-8) {
    dftbMethod_->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
  }
  dftbMethod_->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  dftbMethod_->calculate("TDDFTB reference calculation.");
  tddftbData_ = std::make_unique<TDDFTBData>(dftbMethod_->getTDDFTBData());
}

Core::Calculator& TDDFTBCalculator::getReferenceCalculator() {
  return *dftbMethod_;
}

const Core::Calculator& TDDFTBCalculator::getReferenceCalculator() const {
  return *dftbMethod_;
}

Utils::Settings& TDDFTBCalculator::settings() {
  return *settings_;
}

const Utils::Settings& TDDFTBCalculator::settings() const {
  return *settings_;
}

std::string TDDFTBCalculator::name() const {
  return model;
}

void TDDFTBCalculator::applySettings() {
}

void TDDFTBCalculator::checkMemoryRequirement(int excitationsDim, int numberOfEnergyLevels) {
  double memRequ = 0.;
  double maxMem = settings_->getDouble(Utils::SettingsNames::maxMemory);
  int nAtoms = tddftbData_->AOInfo.getNAtoms();
  int maxSubspaceDim = 0;
  if (numberOfEnergyLevels < 1000) { // Formula for good number of iterations before subspace collapse.
    maxSubspaceDim = static_cast<int>((6 * numberOfEnergyLevels + 50) - pow(numberOfEnergyLevels, 1.2));
    maxSubspaceDim += numberOfEnergyLevels - maxSubspaceDim % numberOfEnergyLevels;
  }
  else {
    maxSubspaceDim = 2 * numberOfEnergyLevels;
  }

  memRequ += static_cast<double>(maxSubspaceDim * excitationsDim * 2);

  for (int i = 0; i < nAtoms; ++i) {
    int nAOsI = tddftbData_->AOInfo.getNOrbitals(i);
    memRequ += static_cast<double>(std::pow(nAOsI, 4));
    for (int j = i + 1; j < nAtoms; ++j) {
      int nAOsJ = tddftbData_->AOInfo.getNOrbitals(j);
      memRequ += static_cast<double>(nAOsI * nAOsI * nAOsJ * nAOsJ);
    }
  }
  memRequ = memRequ * sizeof(double) * 1.e-9;
  getLog().output << "Memory Required: " << memRequ << " GB"
                  << " (specified maximum: " << maxMem << " GB).";
  if (memRequ > maxMem) {
    throw std::runtime_error("This calculation setup would require more memory than specified!"); // TODO setting for
                                                                                                  // allowed memory
  }
}

const Utils::Results& TDDFTBCalculator::calculate() {
  if (!dftbMethod_) {
    throw std::runtime_error("No reference calculator assigned.");
  }
  if (!tddftbData_) {
    throw std::runtime_error("No TDDFTBData was set.");
  }
  if (!dftbMethod_->results().has<Utils::Property::Energy>()) {
    throw InvalidReferenceCalculationException();
  }

  Utils::SpinAdaptedElectronicTransitionResult transitionResult{};

  bool is_unrestricted = tddftbData_->occupation.isUnrestricted();

  std::shared_ptr<Utils::DipoleMatrix> dipoleMatrix;
  if (getReferenceCalculator().results().has<Utils::Property::DipoleMatrixMO>()) {
    dipoleMatrix =
        std::make_shared<Utils::DipoleMatrix>(getReferenceCalculator().results().get<Utils::Property::DipoleMatrixMO>());
  }

  TransitionChargesCalculator tcCalc(tddftbData_->molecularOrbitals, tddftbData_->overlapMatrix, tddftbData_->AOInfo);

  if (!is_unrestricted) { // Here choose if restricted or unrestricted reference
    evalImpl<Utils::Reference::Restricted>(
        transitionResult,
        tcCalc.calculateAtomicTransitionChargeMatrices<Utils::Reference::Restricted>(tddftbData_->occupation), dipoleMatrix);
  }
  else { // If unrestricted
    evalImpl<Utils::Reference::Unrestricted>(
        transitionResult,
        tcCalc.calculateAtomicTransitionChargeMatrices<Utils::Reference::Unrestricted>(tddftbData_->occupation), dipoleMatrix);
  }

  results_.set<Utils::Property::ExcitedStates>(std::move(transitionResult));
  guess_.reset();
  return results_;
}

inline auto orderGuess(std::shared_ptr<TDDFTBCalculator::GuessSpecifier> guess, std::vector<int> orderMap)
    -> std::shared_ptr<TDDFTBCalculator::GuessSpecifier> {
  auto orderedGuess = std::make_shared<TDDFTBCalculator::GuessSpecifier>();
  if (guess->singlet.size() != 0)
    TimeDependentUtils::transformOrder(guess->singlet, orderedGuess->singlet, orderMap, TimeDependentUtils::Direction::To);
  if (guess->triplet.size() != 0)
    TimeDependentUtils::transformOrder(guess->triplet, orderedGuess->triplet, orderMap, TimeDependentUtils::Direction::To);
  if (guess->unrestricted.size() != 0)
    TimeDependentUtils::transformOrder(guess->unrestricted, orderedGuess->unrestricted, orderMap,
                                       TimeDependentUtils::Direction::To);
  return orderedGuess;
}

namespace detail {
auto getLabels(const OrderedInput<Utils::Reference::Restricted>& input) -> std::vector<std::string> {
  return TimeDependentUtils::generateExcitationsLabels(input.excitations());
}
auto getLabels(const OrderedInput<Utils::Reference::Unrestricted>& input) -> std::vector<std::string> {
  return TimeDependentUtils::generateExcitationsLabels(input.excitations(), input.isBeta());
}
} // namespace detail

template<Utils::Reference restrictedness>
void TDDFTBCalculator::evalImpl(Utils::SpinAdaptedElectronicTransitionResult& transitionResult,
                                const Eigen::MatrixXd& transitionCharges, std::shared_ptr<Utils::DipoleMatrix> dipoleMatrix) {
  bool prunedCalculation =
      settings().getString(Utils::SettingsNames::pruneBasis) == Utils::SettingsNames::PruningOptions::energy;

  auto excitations =
      TimeDependentUtils::generateExcitations<restrictedness>(tddftbData_->molecularOrbitals, tddftbData_->occupation);
  Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd> energyDifferenceVector;
  TimeDependentUtils::generateEnergyDifferenceVector(tddftbData_->MOEnergies, tddftbData_->occupation, energyDifferenceVector);
  orderMap_ = TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector);

  OrderedInput<restrictedness> orderedInput(energyDifferenceVector, excitations, transitionCharges, orderMap_);
  if (guess_) {
    guess_ = orderGuess(guess_, orderMap_);
  }

  int numberOfRoots = settings().getInt(Utils::SettingsNames::numberOfEigenstates);
  std::unique_ptr<TDDFTBEigenvalueSolver<restrictedness>> evSolver;

  if (prunedCalculation) {
    auto pruner = BasisPruner<restrictedness>(orderedInput, tddftbData_->gammaMatrix, tddftbData_->spinConstants);
    auto pt = PerturbativeThreshold{settings().getDouble(Utils::SettingsNames::perturbativeThreshold)};
    std::unique_ptr<OrderedInput<restrictedness>> prunedData;
    auto en = EnergyThreshold{settings().getDouble(Utils::SettingsNames::energyThreshold)};

    // Depending on whether the energy threshold is set or not, use energy to
    // infer the number of roots needed or the given number of roots directly.
    if (en.threshold == 0.0) { // Empty settings for threshold has a default 0.0.
      prunedData = std::make_unique<OrderedInput<restrictedness>>(pruner.prune(NumberOfConfigurations{numberOfRoots}, pt));
    }
    else {
      prunedData = std::make_unique<OrderedInput<restrictedness>>(pruner.prune(en, pt));
      numberOfRoots = pruner.getNumberOfRootsUnderThreshold();
    }

    excitations_ = prunedData->excitations();
    transitionResult.transitionLabels = detail::getLabels(*prunedData);

    evSolver = std::make_unique<TDDFTBEigenvalueSolver<restrictedness>>(
        *settings_, tddftbData_->gammaMatrix, tddftbData_->spinConstants, *prunedData, std::move(dipoleMatrix), getLog());

    if (guess_) {
      evSolver->setGuess(pruner.prune(*guess_));
    }
  }
  else { // not pruned
    excitations_ = orderedInput.excitations();
    transitionResult.transitionLabels = detail::getLabels(orderedInput);

    evSolver = std::make_unique<TDDFTBEigenvalueSolver<restrictedness>>(
        *settings_, tddftbData_->gammaMatrix, tddftbData_->spinConstants, orderedInput, std::move(dipoleMatrix), getLog());

    if (guess_) {
      evSolver->setGuess(guess_);
    }
  }

  int initialSubspaceDimension = std::max(settings().getInt(Utils::SettingsNames::initialSubspaceDimension), numberOfRoots);
  solveEigenvalueProblem(transitionResult, evSolver, numberOfRoots, initialSubspaceDimension);
}

template<>
void TDDFTBCalculator::solveEigenvalueProblem<Utils::Reference::Restricted>(
    Utils::SpinAdaptedElectronicTransitionResult& transitionResult,
    const std::unique_ptr<TDDFTBEigenvalueSolver<Utils::Reference::Restricted>>& solver, int numberOfRoots,
    int initialSubspaceDimension) {
  if (isDFTB0(dftbMethod_)) {
    transitionResult.singlet = std::make_shared<Utils::ElectronicTransitionResult>(solver->solveDftb0(numberOfRoots));
    transitionResult.triplet = std::make_shared<Utils::ElectronicTransitionResult>(*transitionResult.singlet);
    transitionResult.triplet->transitionDipoles =
        Eigen::Matrix3Xd::Zero(3, transitionResult.triplet->eigenStates.eigenValues.size());
  }
  else { // SCC
    auto spinBlock = settings_->getString(Utils::SettingsNames::spinBlock);
    if (spinBlock == Utils::SettingsNames::SpinBlocks::singlet || spinBlock == "both") { // If singlet
      transitionResult.singlet = std::make_shared<Utils::ElectronicTransitionResult>(
          solver->solve(numberOfRoots, initialSubspaceDimension, Utils::SpinTransition::Singlet));
    }
    if (spinBlock == Utils::SettingsNames::SpinBlocks::triplet || spinBlock == "both") { // If triplet
      transitionResult.triplet = std::make_shared<Utils::ElectronicTransitionResult>(
          solver->solve(numberOfRoots, initialSubspaceDimension, Utils::SpinTransition::Triplet));
    }
  }
}

template<>
void TDDFTBCalculator::solveEigenvalueProblem<Utils::Reference::Unrestricted>(
    Utils::SpinAdaptedElectronicTransitionResult& transitionResult,
    const std::unique_ptr<TDDFTBEigenvalueSolver<Utils::Reference::Unrestricted>>& solver, int numberOfRoots,
    int initialSubspaceDimension) {
  transitionResult.unrestricted =
      std::make_shared<Utils::ElectronicTransitionResult>(solver->solve(numberOfRoots, initialSubspaceDimension));
}

bool TDDFTBCalculator::isDFTB0(std::shared_ptr<DFTBMethodWrapper> method) const {
  return static_cast<bool>(std::dynamic_pointer_cast<DFTB0MethodWrapper>(method));
}

void TDDFTBCalculator::setGuess(std::shared_ptr<GuessSpecifier> guessVectorMatrix) {
  guess_ = std::move(guessVectorMatrix);
}

// TODO: Now works just for restricted
auto TDDFTBCalculator::getGuess() const -> std::shared_ptr<GuessSpecifier> {
  assert(tddftbData_->occupation.isRestricted());
  if (results_.has<Utils::Property::ExcitedStates>()) {
    GuessSpecifier result;
    const auto& excitedStates = results_.get<Utils::Property::ExcitedStates>();

    int nOcc = tddftbData_->occupation.numberOccupiedRestrictedOrbitals();
    int nVir = tddftbData_->molecularOrbitals.numberOrbitals() - nOcc;

    auto formGuess = [&](std::shared_ptr<Utils::ElectronicTransitionResult> res) -> Eigen::MatrixXd {
      if (res) {
        Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(orderMap_.size(), res->eigenStates.eigenValues.size());
        int index = 0;
        for (const auto& excitation : excitations_) {
          int occ = excitation.occ;
          int vir = excitation.vir - nOcc;
          guess.row(occ * nVir + vir) = res->eigenStates.eigenVectors.row(index++);
        }
        return guess;
      }
      return Eigen::MatrixXd(0, 0);
    };
    result.unrestricted = formGuess(excitedStates.unrestricted);
    result.singlet = formGuess(excitedStates.singlet);
    result.triplet = formGuess(excitedStates.triplet);
    if (result.unrestricted.size() + result.singlet.size() + result.triplet.size() == 0) {
      return std::shared_ptr<GuessSpecifier>();
    }
    else {
      return std::make_shared<GuessSpecifier>(result);
    }
  }
  else {
    throw std::runtime_error("Guess from previous calculation required, but no calculation performed!");
  }
}

TDDFTBCalculator::~TDDFTBCalculator() = default;

} // namespace Sparrow
} // namespace Scine
