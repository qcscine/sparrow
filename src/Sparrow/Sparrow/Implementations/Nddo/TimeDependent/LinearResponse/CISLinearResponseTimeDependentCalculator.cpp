/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal dependencies */
#include "CISLinearResponseTimeDependentCalculator.h"
#include "CISSettings.h"
#include "CISSpinContaminator.h"
#include <Sparrow/Implementations/Exceptions.h>
#include <Sparrow/Implementations/Nddo/NDDOMethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/Global2c2eMatrix.h>
#include <Sparrow/Implementations/TimeDependent/DiagonalPreconditionerEvaluator.h>
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/IO/NativeFilenames.h>
/* External dependencies */
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/Constants.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
#include <Utils/Math/IterativeDiagonalizer/DiagonalizerSettings.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/TimeDependent/TransitionDipoleCalculator.h>
#include <cereal/archives/xml.hpp>
#include <cereal/types/map.hpp>
#include <fstream>

namespace Scine {
namespace Sparrow {
InvalidCalculatorTypeForCIS::InvalidCalculatorTypeForCIS(std::shared_ptr<Scine::Core::Calculator> method) {
  error = "Calculator " + method->name() + " is not an NDDO method and cannot calculate CIS.";
}

void CISLinearResponseTimeDependentCalculator::setReferenceCalculator(std::shared_ptr<Scine::Core::Calculator> method) {
  nddoMethod_ = std::dynamic_pointer_cast<NDDOMethodWrapper>(method);
  if (!nddoMethod_) {
    throw InvalidCalculatorTypeForCIS(method);
  }
  cisData_ = std::make_unique<CISData>(nddoMethod_->getCISData());
}

CISLinearResponseTimeDependentCalculator::CISLinearResponseTimeDependentCalculator() {
  settings_ = std::make_unique<CISSettings>();
}

Utils::Results& CISLinearResponseTimeDependentCalculator::results() {
  return results_;
}

const Utils::Results& CISLinearResponseTimeDependentCalculator::results() const {
  return results_;
}

void CISLinearResponseTimeDependentCalculator::referenceCalculation() {
  if (!nddoMethod_)
    throw MissingReferenceCalculatorException();
  if (nddoMethod_->settings().getDouble(Utils::SettingsNames::selfConsistenceCriterion) > 1e-8) {
    nddoMethod_->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
  }
  nddoMethod_->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  nddoMethod_->calculate("CIS reference calculation.");
  cisData_ = std::make_unique<CISData>(nddoMethod_->getCISData());
}

Core::Calculator& CISLinearResponseTimeDependentCalculator::getReferenceCalculator() {
  return *nddoMethod_;
}

const Core::Calculator& CISLinearResponseTimeDependentCalculator::getReferenceCalculator() const {
  return *nddoMethod_;
}

Utils::Settings& CISLinearResponseTimeDependentCalculator::settings() {
  return *settings_;
}

const Utils::Settings& CISLinearResponseTimeDependentCalculator::settings() const {
  return *settings_;
}

std::string CISLinearResponseTimeDependentCalculator::name() const {
  return model;
}

void CISLinearResponseTimeDependentCalculator::applySettings() {
}

inline Eigen::MatrixXd generateGuess(int nConfigurations, int initialSubspaceDimension,
                                     std::shared_ptr<LinearResponseCalculator::GuessSpecifier> guess,
                                     Utils::SpinTransition spinBlock, Utils::Reference restrictedness) {
  if (!guess) {
    srand(42);
    Eigen::MatrixXd guessVectors = Eigen::MatrixXd::Identity(nConfigurations, initialSubspaceDimension);
    guessVectors.block(0, 0, nConfigurations, initialSubspaceDimension) +=
        0.01 * Eigen::MatrixXd::Random(nConfigurations, initialSubspaceDimension);
    return guessVectors;
  }
  else {
    if (restrictedness == Utils::Reference::Unrestricted) {
      assert(guess->unrestricted.rows() == nConfigurations);
      return guess->unrestricted;
    }
    else if (spinBlock == Utils::SpinTransition::Singlet) {
      assert(guess->singlet.rows() == nConfigurations);
      return guess->singlet;
    }
    else {
      assert(guess->triplet.rows() == nConfigurations);
      return guess->triplet;
    }
  }
}

inline int getNConfigurations(Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> energyDiffVector) {
  return energyDiffVector.restricted.size();
}
inline int getNConfigurations(Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd> energyDiffVector) {
  return energyDiffVector.alpha.size() + energyDiffVector.beta.size();
}

template<Utils::Reference restrictedness>
Utils::ElectronicTransitionResult
CISLinearResponseTimeDependentCalculator::solve(Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd> energyDifferenceVector,
                                                int numberOfEnergyLevels, int initialSubspaceDimension,
                                                Utils::SpinTransition spinBlock) {
  const int nConfigurations = getNConfigurations(energyDifferenceVector);

  // If 0 is given, then give all the energy levels
  if (numberOfEnergyLevels == 0 || numberOfEnergyLevels > nConfigurations) {
    numberOfEnergyLevels = nConfigurations;
  }
  if (initialSubspaceDimension == 0 || initialSubspaceDimension < numberOfEnergyLevels ||
      initialSubspaceDimension > nConfigurations) {
    initialSubspaceDimension = numberOfEnergyLevels;
  }

  checkMemoryRequirement(nConfigurations, numberOfEnergyLevels);

  Utils::NonOrthogonalDavidson diagonalizer(numberOfEnergyLevels, nConfigurations);
  diagonalizer.settings().modifyInt(Utils::initialGuessDimensionOption, initialSubspaceDimension);
  diagonalizer.settings().modifyDouble(Utils::residualNormToleranceOption, settings_->getDouble(convergence));
  diagonalizer.settings().modifyString("gep_algo", settings_->getString("gep_algo"));
  if (settings_->getInt(Utils::SettingsNames::maxDavidsonIterations) != 0) {
    diagonalizer.settings().modifyInt(Utils::SettingsNames::maxDavidsonIterations,
                                      settings_->getInt(Utils::SettingsNames::maxDavidsonIterations));
  }
  diagonalizer.setGuess(generateGuess(nConfigurations, initialSubspaceDimension, guess_, spinBlock, restrictedness));

  diagonalizer.setPreconditionerEvaluator(
      std::make_unique<DiagonalPreconditionerEvaluator>(energyDifferenceVector, OrderTag{}));
  diagonalizer.setSigmaVectorEvaluator(std::make_unique<CISSigmaVectorEvaluator<restrictedness>>(
      *cisData_, excitedStatesParam_, energyDifferenceVector, integralsThresholds_, orderMap_, spinBlock));

  Utils::ElectronicTransitionResult excitedStates;
  excitedStates.eigenStates = diagonalizer.solve(getLog());

  generateTransitionDipoleMoments<restrictedness>(excitedStates, *cisData_, spinBlock);

  //    results.spinContamination = CISSpinContaminator::calculateSpinContaminationOpenShell(
  //        cisData_->molecularOrbitals, results.eigenStates.eigenVectors,
  //        cisData_->occupation.getFilledAlphaOrbitals(), cisData_->occupation.getFilledBetaOrbitals(),
  //        excitationsIndices);

  return excitedStates;
}

void CISLinearResponseTimeDependentCalculator::checkMemoryRequirement(int excitationsDim, int numberOfEnergyLevels) {
  double memRequ = 0.;
  double maxMem = settings_->getDouble(Utils::SettingsNames::maxMemory);
  int nAtoms = cisData_->AOInfo.getNAtoms();
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
    int nAOsI = cisData_->AOInfo.getNOrbitals(i);
    memRequ += static_cast<double>(std::pow(nAOsI, 4));
    for (int j = i + 1; j < nAtoms; ++j) {
      int nAOsJ = cisData_->AOInfo.getNOrbitals(j);
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

const Utils::Results& CISLinearResponseTimeDependentCalculator::calculate() {
  prepareIntegralScreening();
  if (!nddoMethod_) {
    throw std::runtime_error("No reference calculator assigned.");
  }
  if (!cisData_) {
    throw std::runtime_error("No CISData was set.");
  }
  if (!nddoMethod_->results().has<Utils::Property::Energy>()) {
    throw InvalidReferenceCalculationException();
  }

  Utils::SpinAdaptedElectronicTransitionResult transitionResult{};

  bool is_unrestricted = cisData_->occupation.isUnrestricted();
  int numberOfRoots = settings().getInt(Utils::SettingsNames::numberOfEigenstates);
  int initialSubspaceDimension = settings().getInt(Utils::SettingsNames::initialSubspaceDimension);

  if (!is_unrestricted) { // Here choose if restricted or unrestricted reference
    Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> energyDifferenceVector;
    TimeDependentUtils::generateEnergyDifferenceVector(cisData_->MOEnergies, cisData_->occupation, energyDifferenceVector);
    orderMap_ = TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector);

    TimeDependentUtils::transformOrder(
        TimeDependentUtils::generateExcitationsLabels(TimeDependentUtils::generateExcitations<Utils::Reference::Restricted>(
            cisData_->molecularOrbitals, cisData_->occupation)),
        transitionResult.transitionLabels, orderMap_, TimeDependentUtils::Direction::To);

    auto spinBlock = settings_->getString(Utils::SettingsNames::spinBlock);

    if (spinBlock == Utils::SettingsNames::SpinBlocks::singlet || spinBlock == "both") { // If singlet
      setExcitedStatesParam(Utils::Reference::Restricted, Utils::SpinTransition::Singlet);
      transitionResult.singlet = std::make_unique<Utils::ElectronicTransitionResult>(
          solve(energyDifferenceVector, numberOfRoots, initialSubspaceDimension, Utils::SpinTransition::Singlet));
    }
    if (spinBlock == Utils::SettingsNames::SpinBlocks::triplet || spinBlock == "both") { // If triplet
      setExcitedStatesParam(Utils::Reference::Restricted, Utils::SpinTransition::Triplet);
      transitionResult.triplet = std::make_unique<Utils::ElectronicTransitionResult>(
          solve(energyDifferenceVector, numberOfRoots, initialSubspaceDimension, Utils::SpinTransition::Triplet));
    }
  }
  else { // If unrestricted
    Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd> energyDifferenceVector;
    TimeDependentUtils::generateEnergyDifferenceVector(cisData_->MOEnergies, cisData_->occupation, energyDifferenceVector);
    orderMap_ = TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector);
    TimeDependentUtils::transformOrder(TimeDependentUtils::generateExcitationsLabels(
                                           TimeDependentUtils::generateExcitations<Utils::Reference::Unrestricted>(
                                               cisData_->molecularOrbitals, cisData_->occupation)),
                                       transitionResult.transitionLabels, orderMap_, TimeDependentUtils::Direction::To);
    setExcitedStatesParam(Utils::Reference::Unrestricted, {});
    transitionResult.unrestricted = std::make_unique<Utils::ElectronicTransitionResult>(
        solve(energyDifferenceVector, numberOfRoots, initialSubspaceDimension));
  }

  results_.set<Utils::Property::ExcitedStates>(std::move(transitionResult));
  guess_.reset();
  return results_;
}

template<Utils::Reference restrictedness>
void CISLinearResponseTimeDependentCalculator::generateTransitionDipoleMoments(Utils::ElectronicTransitionResult& excitedStatesResults,
                                                                               const CISData& cisData,
                                                                               Utils::SpinTransition spinBlock) const {
  assert(!orderMap_.empty());
  // Triplets have a transition dipole moment of 0, the function transitionDipoleMoment just multiplies by 2 in the
  // restricted case so care must be taken.
  if (getReferenceCalculator().results().has<Utils::Property::DipoleMatrixMO>() && spinBlock == Utils::SpinTransition::Singlet) {
    std::vector<Utils::Excitation> orderedExcitation;
    TimeDependentUtils::transformOrder(TimeDependentUtils::flatten(TimeDependentUtils::generateExcitations<restrictedness>(
                                           cisData.molecularOrbitals, cisData.occupation)),
                                       orderedExcitation, orderMap_, TimeDependentUtils::Direction::To);

    excitedStatesResults.transitionDipoles = Utils::TransitionDipoleCalculator::calculate<restrictedness>(
        getReferenceCalculator().results().get<Utils::Property::DipoleMatrixMO>(),
        excitedStatesResults.eigenStates.eigenVectors, orderedExcitation);
  }
  else {
    excitedStatesResults.transitionDipoles = Eigen::Matrix3Xd::Zero(3, excitedStatesResults.eigenStates.eigenValues.size());
  }
}

void CISLinearResponseTimeDependentCalculator::prepareIntegralScreening() {
  int nAtoms = cisData_->AOInfo.getNAtoms();
  integralsThresholds_.clear();
  integralsThresholds_.resize(nAtoms);
  for (int atomI = 0; atomI < nAtoms; ++atomI) {
    for (int atomJ = atomI + 1; atomJ < nAtoms; ++atomJ) {
      const auto& twoCenterBlock = cisData_->twoCenterIntegrals.get(atomI, atomJ)->getGlobalMatrix();
      double maxIntegral = std::max(std::abs(twoCenterBlock.maxCoeff()), std::abs(twoCenterBlock.minCoeff()));
      integralsThresholds_[atomI].insert({maxIntegral, atomJ});
    }
  }
}

void CISLinearResponseTimeDependentCalculator::setExcitedStatesParam(Utils::Reference restrictedness,
                                                                     Utils::SpinTransition spinBlock) {
  std::string path = settings_->getString(Utils::SettingsNames::excitedStatesParamFile);
  std::map<std::string, ExcitedStatesParam> readExcParam = {};

  if (!path.empty()) {
    std::ifstream fs(path);
    if (!fs.is_open())
      throw Utils::Methods::ParameterFileCannotBeOpenedException(path);
    cereal::XMLInputArchive archive(fs);
    archive(readExcParam);
    if (restrictedness == Utils::Reference::Restricted) {
      if (spinBlock == Utils::SpinTransition::Singlet) {
        excitedStatesParam_ = readExcParam["S"];
      }
      if (spinBlock == Utils::SpinTransition::Triplet) {
        excitedStatesParam_ = readExcParam["T"];
      }
    }
    else if (restrictedness == Utils::Reference::Unrestricted) {
      excitedStatesParam_ = readExcParam["U"];
    }
    else
      throw std::runtime_error("CISLinearResponseTimeDependentCalulator: Invalid Reference");
  }
  else {
    excitedStatesParam_ = {1., 1., 1.};
  }
}

void CISLinearResponseTimeDependentCalculator::setGuess(std::shared_ptr<GuessSpecifier> guessVectorMatrix) {
  guess_ = std::move(guessVectorMatrix);
}

auto CISLinearResponseTimeDependentCalculator::getGuess() const -> std::shared_ptr<GuessSpecifier> {
  if (results_.has<Utils::Property::ExcitedStates>()) {
    auto guess = std::make_shared<GuessSpecifier>();
    const auto& excitedStates = results_.get<Utils::Property::ExcitedStates>();
    if (excitedStates.singlet) {
      guess->singlet = excitedStates.singlet->eigenStates.eigenVectors;
    }
    if (excitedStates.triplet) {
      guess->triplet = excitedStates.triplet->eigenStates.eigenVectors;
    }
    if (excitedStates.singlet) {
      guess->unrestricted = excitedStates.unrestricted->eigenStates.eigenVectors;
    }
    return guess;
  }
  else {
    return std::shared_ptr<GuessSpecifier>();
  }
}

} // namespace Sparrow
} // namespace Scine
