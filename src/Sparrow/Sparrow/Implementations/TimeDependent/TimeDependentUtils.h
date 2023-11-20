/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_TIMEDEPENDENTUTILS_H
#define SPARROW_TIMEDEPENDENTUTILS_H

#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/TimeDependent/TransitionDipoleCalculator.h>
#include <Eigen/Core>
#include <cereal/cereal.hpp>
#include <map>

namespace Scine {

namespace Sparrow {

struct ExcitedStatesParam {
  double c1 = 1, c2 = 1, d = 1;

  template<class Archive>
  void serialize(Archive& archive) {
    archive(CEREAL_NVP(c1), CEREAL_NVP(c2), CEREAL_NVP(d));
  }
};

/**
 * @brief Class containing utility functions for excites-states calculations such as (NDDO-)CIS or TD-DFT(B).
 * @class TimeDependentUtils @file TimeDependentUtils.h
 * This class contains static functions to calculate occupied-virtual energy difference vectors,
 * generate excitation labels (occ -> vir) for restricted and unrestricted reference calculations.
 */
class TimeDependentUtils {
 public:
  /**
   * @brief Calculates the occupied->virtual energy difference vector.
   * @pre This function assumes an electronic occupation according to the Aufbau principle
   * (no fractional occupation).
   * @param nOccupiedOrbitals The number of occupied orbitals.
   * @param nVirtualOrbitals The number of virtual orbitals.
   * @param energies The list of molecular orbital energies.
   */
  static Eigen::VectorXd generateEnergyDifferenceVector(int nOccupiedOrbitals, int nVirtualOrbitals,
                                                        const Utils::SingleParticleEnergies& energies);
  /**
   * @brief Function generating a vector of "occ->vir" excitation labels.
   * @pre This function assumes an occupation according to the Aufbau principle
   * (no fractional occupation).
   * @param excitations A SpinAdaptedContainer vector containing the occupied->virtual (alpha/beta for unrestricted)
   * excitations in standard order (vir fast index, occ slow index).
   * @tparam restrictedness Flags whether it is a Restricted or Unrestricted reference calculation.
   * @return A vector with the excitations in same order as excitations, alpha/beta concatenated.
   */
  template<Utils::Reference restrictedness>
  static std::vector<std::string>
  generateExcitationsLabels(const Utils::SpinAdaptedContainer<restrictedness, std::vector<Utils::Excitation>>& excitations);

  /**
   * @brief Function generating a vector of "occ->vir" excitation labels.
   * @pre This function assumes an occupation according to the Aufbau principle
   * (no fractional occupation).
   * @param excitations A vector containing the occupied->virtual (alpha/beta for unrestricted) excitations in arbitrary
   *        order.
   * @param isBeta A bool vector indicating whether the excitation is alpha or beta.
   * @return A vector with the excitations in same order as excitations.
   */
  static std::vector<std::string> generateExcitationsLabels(const std::vector<Utils::Excitation>& excitations,
                                                            const Eigen::Matrix<bool, -1, 1>& isBeta = {});

  /**
   * @brief Calculated the occupied->virtual energy difference vector.
   * @param mos The molecular orbitals, needed to know the number of occupied and virtual orbitals.
   * @param energies The list of molecular orbital energies.
   * @param occupation The electronic, usually according to Aufbau principle. Here only integer occupations supported.
   * @tparam restrictedness Flags whether it is a Restricted or Unrestricted reference calculation.
   */
  template<Utils::Reference restrictedness>
  static void generateEnergyDifferenceVector(const Utils::SingleParticleEnergies& energies,
                                             const Utils::LcaoUtils::ElectronicOccupation& occupation,
                                             Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd>& energyDifferenceVector);

  /**
   * @brief Fills a vector of virtual orbital indices having just size and occupied indices.
   * @param occupiedOrbitalsIndices A vector<int> containing the indices of the occupied orbitals.
   * @param virtualOrbitalsIndices A vector<int> with size equal to the number of virtual orbitals,
   *                               given by reference as its elements will be filled with the
   *                               virtual orbital indices.
   */
  static void generateVirtualOrbitalIndices(const std::vector<int>& occupiedOrbitalsIndices,
                                            std::vector<int>& virtualOrbitalsIndices);
  /**
   * @brief Given occupied and virtual orbital indices, generates a vector of occ->vir excitations.
   * @param occupiedOrbitalsIndices A vector containing the index corresponding to occupied orbitals.
   * @param virtualOrbitalsIndices A vector containing the index corresponding to virtual orbitals.
   * @return A vector of excitations, objects with an occupied and a virtual index.
   */
  static std::vector<Utils::Excitation> generateExcitations(const std::vector<int>& occupiedOrbitalsIndices,
                                                            const std::vector<int>& virtualOrbitalsIndices);
  /**
   * @brief Given occupied and virtual orbital indices, generates a vector of occ->vir excitations.
   * @param occupiedOrbitalsIndices A vector containing the index corresponding to occupied orbitals.
   * @param virtualOrbitalsIndices A vector containing the index corresponding to virtual orbitals.
   * @tparam restrictedness Flags whether it is a Restricted or Unrestricted reference calculation.
   * @return A SpinAdaptedContained of excitations, objects with an occupied and a virtual index.
   *         If restricted, just contains one, if unrestricted, then contains alpha and beta.
   */
  template<Utils::Reference restrictedness>
  static Utils::SpinAdaptedContainer<restrictedness, std::vector<Utils::Excitation>>
  generateExcitations(const Utils::MolecularOrbitals& mos, const Utils::LcaoUtils::ElectronicOccupation& occupation);

  /**
   * @brief Generates a mapping from the "standard" excitation order to increasing energy order.
   * @return A vector with the indices in the new order.
   * For example, let's assume there is an orbital energy difference vector en with elements:
   * en[0] = 3.5
   * en[1] = 2.2
   * en[2] = 6.4
   * en[3] = 1.0
   *
   * then the resulting vector res will contain
   * res[0] = 3
   * res[1] = 1
   * res[2] = 0
   * res[3] = 2
   *
   * So, to map an arbitrary vector standardV FROM the "standard" order TO the
   * increasing energy ordering in a resut vector increasingEnergyV:
   *
   * for (int i = 0; i < res.size(); ++i)
   *   increasingEnergyV[i] = standardV[res[i]];
   *
   * And, to map an arbitrary vector increasingEnergyV FROM the increasing energy order TO the
   * "standard" ordering in a resut vector standardV:
   *
   * for (int i = 0; i < res.size(); ++i)
   *   standardV[res[i]] = increasingEnergyV[i];
   */
  template<Utils::Reference restrictedness>
  static std::vector<int>
  generateEnergyOrderMap(const Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd>& energyDifferenceVector);

  enum class Direction { To, From };
  /**
   * @brief Transforms a vector or matrix from one ordering to another, or back.
   * @param toTransform The vector or matrix whose rows to transform.
   * @param orderVector The vector with the ordering.
   *                    For example, let's assume there is an orbital energy
   *                    difference vector en with elements:
   *                    en[0] = 3.5
   *                    en[1] = 2.2
   *                    en[2] = 6.4
   *                    en[3] = 1.0
   *
   *                    then the resulting mapping vector res will contain
   *                    res[0] = 3
   *                    res[1] = 1
   *                    res[2] = 0
   *                    res[3] = 2
   *                    This is also the vector returned by
   *                    TimeDependentUtils::generateEnergyOrderMap()
   * @param direction The direction of the transformation.
   *              - Direction::To :
   *                transform a vector as implied by orderVector, i.e.
   *                in an the ordering according to energy difference, from
   *                the "standard" order to the order of increasing energy.
   *              - Direction::From :
   *                The opposite direction of transformation.
   * @return A vector transformed into the new order defined in the argument.
   *
   * So, to map an arbitrary vector standardV FROM the "standard" order TO the
   * increasing energy ordering in a resut vector increasingEnergyV:
   *
   * Eigen::VectorXd increasingEnergyV = transformOrder(standardV, orderVec, Direction::To);
   * and to transform it back:
   * Eigen::VectorXd standardV2 = transformOrder(increasingEnergyV, orderVec, Direction::From);
   * Then standardV and standardV2 should be equal.
   *
   */
  template<typename Derived, typename DerivedResult>
  static void transformOrder(const Eigen::MatrixBase<Derived>& toTransform, const Eigen::MatrixBase<DerivedResult>& result,
                             const std::vector<int>& orderVector, Direction direction);
  /**
   * @brief Transforms a std::vector from one ordering to another, or back.
   */
  template<typename Type>
  static void transformOrder(const std::vector<Type>& toTransform, std::vector<Type>& result,
                             const std::vector<int>& orderVector, Direction direction);

  /**
   * @brief Transforms the SpinAdaptedContainer to a Eigen::VectorXd by extracting its "restricted" member.
   */
  static auto flatten(const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& in) -> Eigen::VectorXd;
  /**
   * @brief Transforms the SpinAdaptedContainer to a Eigen::VectorXd by concatenating its "alpha" and "beta" member.
   */
  static auto flatten(const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& in) -> Eigen::VectorXd;

  /**
   * @brief Transforms the SpinAdaptedContainer to a std::vector<T> by extracting its "restricted" member.
   */
  template<typename T>
  static auto flatten(Utils::SpinAdaptedContainer<Utils::Reference::Restricted, std::vector<T>> in) -> std::vector<T>;

  /**
   * @brief Transforms the SpinAdaptedContainer to a std::vector<T> by concatenating its "alpha" and "beta" member.
   */
  template<typename T>
  static auto flatten(Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<T>> in) -> std::vector<T>;
};

inline auto TimeDependentUtils::flatten(const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& in)
    -> Eigen::VectorXd {
  return in.restricted;
}
inline auto TimeDependentUtils::flatten(const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& in)
    -> Eigen::VectorXd {
  Eigen::VectorXd out(in.alpha.size() + in.beta.size());
  out << in.alpha, in.beta;
  return out;
}

template<typename T>
inline auto TimeDependentUtils::flatten(Utils::SpinAdaptedContainer<Utils::Reference::Restricted, std::vector<T>> in)
    -> std::vector<T> {
  return in.restricted;
}
template<typename T>
inline auto TimeDependentUtils::flatten(Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<T>> in)
    -> std::vector<T> {
  std::vector<T> out;
  out.reserve(in.alpha.size() + in.beta.size());
  std::copy(in.alpha.begin(), in.alpha.end(), std::back_inserter(out));
  std::copy(in.beta.begin(), in.beta.end(), std::back_inserter(out));
  return out;
}

inline Eigen::VectorXd TimeDependentUtils::generateEnergyDifferenceVector(int nOccupiedOrbitals, int nVirtualOrbitals,
                                                                          const Utils::SingleParticleEnergies& energies) {
  auto const& energyLevels = energies.getRestrictedEnergies();

  Eigen::VectorXd energyDifferenceVector(nOccupiedOrbitals * nVirtualOrbitals);
  for (int occ = 0; occ < nOccupiedOrbitals; ++occ) {
    for (int vir = 0; vir < nVirtualOrbitals; ++vir) {
      energyDifferenceVector(nVirtualOrbitals * occ + vir) = energyLevels[nOccupiedOrbitals + vir] - energyLevels[occ];
    }
  }
  return energyDifferenceVector;
}

template<>
inline void TimeDependentUtils::generateEnergyDifferenceVector<Utils::Reference::Unrestricted>(
    const Utils::SingleParticleEnergies& energies, const Utils::LcaoUtils::ElectronicOccupation& occupation,
    Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& energyDifferenceVector) {
  const auto& alphaEnergies = energies.getAlphaEnergies();
  const auto& betaEnergies = energies.getBetaEnergies();
  auto const nAlphaMOs = alphaEnergies.size();
  auto const nBetaMOs = betaEnergies.size();
  const std::vector<int>& occupiedAlphaOrbitals = occupation.getFilledAlphaOrbitals();
  const std::vector<int>& occupiedBetaOrbitals = occupation.getFilledBetaOrbitals();
  auto const nOccAlpha = occupiedAlphaOrbitals.size();
  auto const nOccBeta = occupiedBetaOrbitals.size();
  auto const nVirAlpha = nAlphaMOs - nOccAlpha;
  auto const nVirBeta = nBetaMOs - nOccBeta;

  energyDifferenceVector.alpha.resize(nOccAlpha * nVirAlpha);
  energyDifferenceVector.beta.resize(nOccBeta * nVirBeta);

  int alphaIter = 0;
  for (auto const& occ : occupiedAlphaOrbitals) {
    for (unsigned int mo = 0; mo < nAlphaMOs; ++mo) {
      if ((std::find(occupiedAlphaOrbitals.begin(), occupiedAlphaOrbitals.end(), mo) == occupiedAlphaOrbitals.end())) {
        energyDifferenceVector.alpha[alphaIter] = alphaEnergies[mo] - alphaEnergies[occ];
        alphaIter++;
      }
    }
  }
  int betaIter = 0;
  for (auto const& occ : occupiedBetaOrbitals) {
    for (unsigned int mo = 0; mo < nBetaMOs; ++mo) {
      if ((std::find(occupiedBetaOrbitals.begin(), occupiedBetaOrbitals.end(), mo) == occupiedBetaOrbitals.end())) {
        energyDifferenceVector.beta[betaIter] = betaEnergies[mo] - betaEnergies[occ];
        betaIter++;
      }
    }
  }
}

template<>
inline void TimeDependentUtils::generateEnergyDifferenceVector<Utils::Reference::Restricted>(
    const Utils::SingleParticleEnergies& energies, const Utils::LcaoUtils::ElectronicOccupation& occupation,
    Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& energyDifferenceVector) {
  std::vector<int> occupiedOrbitals = occupation.getFilledRestrictedOrbitals();
  const auto& energiesLevels = energies.getRestrictedEnergies();
  auto const nMOs = energiesLevels.size();
  auto const nOcc = occupiedOrbitals.size();
  auto const nVir = nMOs - nOcc;

  energyDifferenceVector.restricted.resize(nOcc * nVir);

  int iter = 0;
  for (int occ : occupiedOrbitals) {
    for (unsigned int mo = 0; mo < nMOs; mo++) {
      if ((std::find(occupiedOrbitals.begin(), occupiedOrbitals.end(), mo) == occupiedOrbitals.end())) {
        energyDifferenceVector.restricted[iter] = energiesLevels[mo] - energiesLevels[occ];
        iter++;
      }
    }
  }
}

inline void TimeDependentUtils::generateVirtualOrbitalIndices(const std::vector<int>& occupiedOrbitalsIndices,
                                                              std::vector<int>& virtualOrbitalsIndices) {
  unsigned int nMOs = occupiedOrbitalsIndices.size() + virtualOrbitalsIndices.size();
  unsigned int iterEmpty = 0;
  unsigned int iterFilled = 0;
  for (int i = 0; i < static_cast<int>(nMOs); ++i) {
    if (iterFilled < occupiedOrbitalsIndices.size() && i == occupiedOrbitalsIndices[iterFilled]) {
      ++iterFilled;
    }
    else {
      virtualOrbitalsIndices[iterEmpty] = i;
      ++iterEmpty;
    }
  }
}

inline std::vector<Utils::Excitation> TimeDependentUtils::generateExcitations(const std::vector<int>& occupiedOrbitalsIndices,
                                                                              const std::vector<int>& virtualOrbitalsIndices) {
  std::vector<Utils::Excitation> excitations(occupiedOrbitalsIndices.size() * virtualOrbitalsIndices.size());

  int iter = 0;
  for (int occ : occupiedOrbitalsIndices) {
    for (int vir : virtualOrbitalsIndices) {
      excitations[iter].occ = occ;
      excitations[iter].vir = vir;
      iter++;
    }
  }
  return excitations;
}
template<>
inline Utils::SpinAdaptedContainer<Utils::Reference::Restricted, std::vector<Utils::Excitation>>
TimeDependentUtils::generateExcitations<Utils::Reference::Restricted>(const Utils::MolecularOrbitals& mos,
                                                                      const Utils::LcaoUtils::ElectronicOccupation& occupation) {
  int nMOs = mos.restrictedMatrix().cols();
  const auto& occupiedOrbitals = occupation.getFilledRestrictedOrbitals();
  std::vector<int> virtualOrbitals(nMOs - occupiedOrbitals.size());

  generateVirtualOrbitalIndices(occupiedOrbitals, virtualOrbitals);
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, std::vector<Utils::Excitation>> excitations;
  excitations.restricted.resize(occupiedOrbitals.size() * virtualOrbitals.size());

  int iter = 0;
  for (int occ : occupiedOrbitals) {
    for (int vir : virtualOrbitals) {
      excitations.restricted[iter].occ = occ;
      excitations.restricted[iter].vir = vir;
      iter++;
    }
  }
  return excitations;
}

template<>
inline Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<Utils::Excitation>>
TimeDependentUtils::generateExcitations<Utils::Reference::Unrestricted>(const Utils::MolecularOrbitals& mos,
                                                                        const Utils::LcaoUtils::ElectronicOccupation& occupation) {
  int nMOs = mos.alphaMatrix().cols();
  const auto& occupiedAlphaOrbitals = occupation.getFilledAlphaOrbitals();
  const auto& occupiedBetaOrbitals = occupation.getFilledBetaOrbitals();
  std::vector<int> virtualAlphaOrbitals(nMOs - occupiedAlphaOrbitals.size());
  std::vector<int> virtualBetaOrbitals(nMOs - occupiedBetaOrbitals.size());

  generateVirtualOrbitalIndices(occupiedAlphaOrbitals, virtualAlphaOrbitals);
  generateVirtualOrbitalIndices(occupiedBetaOrbitals, virtualBetaOrbitals);
  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<Utils::Excitation>> excitations;
  excitations.alpha.resize(occupiedAlphaOrbitals.size() * virtualAlphaOrbitals.size());
  excitations.beta.resize(occupiedBetaOrbitals.size() * virtualBetaOrbitals.size());

  int iter = 0;
  for (int occ : occupiedAlphaOrbitals) {
    for (int vir : virtualAlphaOrbitals) {
      excitations.alpha[iter].occ = occ;
      excitations.alpha[iter].vir = vir;
      iter++;
    }
  }
  iter = 0;
  for (int occ : occupiedBetaOrbitals) {
    for (int vir : virtualBetaOrbitals) {
      excitations.beta[iter].occ = occ;
      excitations.beta[iter].vir = vir;
      iter++;
    }
  }
  return excitations;
}

namespace detail {
struct BetaExcitation {
  bool value = true;
};

inline auto occToVirLabel(int occupiedIndex, int virtualIndex) -> std::string {
  return std::to_string(occupiedIndex) + " -> " + std::to_string(virtualIndex);
}

inline auto occToVirLabel(int occupiedIndex, int virtualIndex, BetaExcitation isBeta) -> std::string {
  std::string spinLabel = isBeta.value ? "b" : "a";
  return std::to_string(occupiedIndex) + spinLabel + " -> " + std::to_string(virtualIndex) + spinLabel;
}
} // namespace detail

template<>
inline std::vector<std::string> TimeDependentUtils::generateExcitationsLabels<Utils::Reference::Restricted>(
    const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, std::vector<Utils::Excitation>>& excitations) {
  int nExc = excitations.restricted.size();
  std::vector<std::string> labels;

  labels.reserve(nExc);

  for (int i = 0; i < nExc; ++i) {
    int indexOccupied = excitations.restricted[i].occ;
    int indexVirtual = excitations.restricted[i].vir;
    labels.push_back(detail::occToVirLabel(indexOccupied, indexVirtual));
  }
  return labels;
}

template<>
inline std::vector<std::string> TimeDependentUtils::generateExcitationsLabels<Utils::Reference::Unrestricted>(
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<Utils::Excitation>>& excitations) {
  int nAlphaExc = excitations.alpha.size();
  int nBetaExc = excitations.beta.size();

  std::vector<std::string> labels;
  labels.reserve(nAlphaExc + nBetaExc);

  for (int i = 0; i < nAlphaExc; ++i) {
    int indexOccupied = excitations.alpha[i].occ;
    int indexVirtual = excitations.alpha[i].vir;
    labels.push_back(detail::occToVirLabel(indexOccupied, indexVirtual, detail::BetaExcitation{false}));
  }

  for (int i = 0; i < nBetaExc; ++i) {
    int indexOccupied = excitations.beta[i].occ;
    int indexVirtual = excitations.beta[i].vir;
    labels.push_back(detail::occToVirLabel(indexOccupied, indexVirtual, detail::BetaExcitation{true}));
  }
  return labels;
}

inline std::vector<std::string> TimeDependentUtils::generateExcitationsLabels(const std::vector<Utils::Excitation>& excitations,
                                                                              const Eigen::Matrix<bool, -1, 1>& isBeta) {
  assert(static_cast<Eigen::Index>(excitations.size()) == isBeta.size() || isBeta.size() == 0);
  std::vector<std::string> labels;
  labels.reserve(excitations.size());

  for (unsigned int i = 0; i < excitations.size(); ++i) {
    int occ = excitations[i].occ;
    int vir = excitations[i].vir;
    if (isBeta.size() != 0) {
      labels.push_back(detail::occToVirLabel(occ, vir, detail::BetaExcitation{isBeta(i)}));
    }
    else {
      labels.push_back(detail::occToVirLabel(occ, vir));
    }
  }
  return labels;
}

template<>
inline std::vector<int> TimeDependentUtils::generateEnergyOrderMap<Utils::Reference::Restricted>(
    const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& energyDifferenceVector) {
  std::multimap<double, int> energyOrderMap;
  std::vector<int> orderedIndices;
  const Eigen::VectorXd& energies = energyDifferenceVector.restricted;
  orderedIndices.reserve(energies.size());
  for (int i = 0; i < energies.size(); ++i) {
    energyOrderMap.insert({energies(i), i});
  }
  for (const auto& element : energyOrderMap) {
    orderedIndices.push_back(element.second);
  }
  return orderedIndices;
}

template<>
inline std::vector<int> TimeDependentUtils::generateEnergyOrderMap<Utils::Reference::Unrestricted>(
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& energyDifferenceVector) {
  std::multimap<double, int> energyOrderMap;
  std::vector<int> orderedIndices;
  const Eigen::VectorXd& alpha = energyDifferenceVector.alpha;
  const Eigen::VectorXd& beta = energyDifferenceVector.beta;
  Eigen::VectorXd concatenated(alpha.size() + beta.size());
  concatenated << alpha, beta;
  orderedIndices.reserve(concatenated.size());
  for (int i = 0; i < concatenated.size(); ++i) {
    energyOrderMap.insert({concatenated(i), i});
  }
  for (const auto& element : energyOrderMap) {
    orderedIndices.push_back(element.second);
  }
  return orderedIndices;
}

template<typename Derived, typename DerivedResult>
inline void TimeDependentUtils::transformOrder(const Eigen::MatrixBase<Derived>& toTransform,
                                               const Eigen::MatrixBase<DerivedResult>& result,
                                               const std::vector<int>& orderVector, Direction direction) {
  assert(direction == Direction::To || direction == Direction::From);
  assert(static_cast<Eigen::Index>(orderVector.size()) == toTransform.rows());
  auto& castResult = const_cast<Eigen::MatrixBase<DerivedResult>&>(result).derived();
  castResult.resize(toTransform.rows(), toTransform.cols());
  if (direction == Direction::To) {
    for (unsigned int i = 0; i < orderVector.size(); ++i) {
      castResult.row(i) = toTransform.row(orderVector[i]);
    }
  }
  else { // Direction::From
    for (unsigned int i = 0; i < orderVector.size(); ++i) {
      castResult.row(orderVector[i]) = toTransform.row(i);
    }
  }
}

template<typename Type>
inline void TimeDependentUtils::transformOrder(const std::vector<Type>& toTransform, std::vector<Type>& result,
                                               const std::vector<int>& orderVector, Direction direction) {
  assert(direction == Direction::To || direction == Direction::From);
  assert(static_cast<Eigen::Index>(orderVector.size()) == static_cast<Eigen::Index>(toTransform.size()));
  result.resize(toTransform.size());
  if (direction == Direction::To) {
    for (unsigned int i = 0; i < orderVector.size(); ++i) {
      result[i] = toTransform[orderVector[i]];
    }
  }
  else { // Direction::From
    for (unsigned int i = 0; i < orderVector.size(); ++i) {
      result[orderVector[i]] = toTransform[i];
    }
  }
}

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_TIMEDEPENDENTUTILS_H
