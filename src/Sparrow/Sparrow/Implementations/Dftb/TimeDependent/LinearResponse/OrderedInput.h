/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ORDEREDINPUT_H
#define SPARROW_ORDEREDINPUT_H

#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>

namespace Scine {
namespace Sparrow {
namespace detail {
/**
 * @brief Result of the pruning procedure. Restricted specialization.
 * Note: all quantities are sorted in increasing energy order.
 */
struct RestrictedBase {
 protected:
  RestrictedBase() = default;
  RestrictedBase(const std::vector<int>& /*orderMap*/) {
  }
};
/**
 * @brief Result of the pruning procedure. Unrestricted specialization.
 * Has an additional "isBeta" bool vector.
 * Note: all quantities are sorted in increasing energy order.
 */
struct UnrestrictedBase {
  auto isBeta() const -> const Eigen::Matrix<bool, -1, 1>& {
    return isBeta_;
  }
  auto isBeta() -> Eigen::Matrix<bool, -1, 1>& {
    return isBeta_;
  }

 protected:
  UnrestrictedBase() = default;

  UnrestrictedBase(const std::vector<int>& orderMap) {
    assert(!orderMap.empty());
    int nTransitions = orderMap.size();
    Eigen::Matrix<bool, -1, 1> unorderedIsBeta = Eigen::Matrix<bool, -1, 1>::Constant(nTransitions, false);
    // Same number of alpha and beta transitions
    unorderedIsBeta.tail(nTransitions / 2).array() = true;

    TimeDependentUtils::transformOrder(unorderedIsBeta, isBeta_, orderMap, TimeDependentUtils::Direction::To);
  }

 private:
  Eigen::Matrix<bool, Eigen::Dynamic, 1> isBeta_;
};
} // namespace detail

/**
 * @brief Ordered input for the TDDFTB eigenvalue solver. Contains quantities in increasing energetic order.
 * @class OrderedInput
 * @tparam restrictedness Decides whether a restricted or unrestricted reference calculation was used.
 * Note: all quantities are sorted in increasing energy order.
 * Note on conditional inheritance:
 * std::conditional<SomeBool, A, B>::type gives A if SomeBool is true, B otherwise.
 * So, the class derives from the (empty) detail::RestrictedBase struct if restrictedness
 * is Restricted, and from the struct providing the isBeta member if restrictedness is
 * Unrestricted.
 */
template<Utils::Reference restrictedness>
class OrderedInput
  : public std::conditional<restrictedness == Utils::Reference::Restricted, detail::RestrictedBase, detail::UnrestrictedBase>::type {
 public:
  using Base =
      typename std::conditional<restrictedness == Utils::Reference::Restricted, detail::RestrictedBase, detail::UnrestrictedBase>::type;
  OrderedInput() = default;
  OrderedInput(const Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd>& energyDifferenceVector,
               const Utils::SpinAdaptedContainer<restrictedness, std::vector<Utils::Excitation>>& excitations,
               const Eigen::MatrixXd& transitionCharges, const std::vector<int>& orderMap)
    : Base(orderMap) {
    TimeDependentUtils::transformOrder(TimeDependentUtils::flatten(energyDifferenceVector), energyDifferences_,
                                       orderMap, TimeDependentUtils::Direction::To);
    TimeDependentUtils::transformOrder(TimeDependentUtils::flatten(excitations), excitations_, orderMap,
                                       TimeDependentUtils::Direction::To);
    TimeDependentUtils::transformOrder(transitionCharges, transitionCharges_, orderMap, TimeDependentUtils::Direction::To);
  }

  auto energyDifferences() const -> const Eigen::VectorXd& {
    return energyDifferences_;
  }
  auto energyDifferences() -> Eigen::VectorXd& {
    return energyDifferences_;
  }
  auto excitations() const -> const std::vector<Utils::Excitation>& {
    return excitations_;
  }
  auto excitations() -> std::vector<Utils::Excitation>& {
    return excitations_;
  }
  auto transitionCharges() const -> const Eigen::MatrixXd& {
    return transitionCharges_;
  }
  auto transitionCharges() -> Eigen::MatrixXd& {
    return transitionCharges_;
  }

 private:
  Eigen::VectorXd energyDifferences_;
  std::vector<Utils::Excitation> excitations_;
  Eigen::MatrixXd transitionCharges_;
};
} // namespace Sparrow
} // namespace Scine

#endif //  SPARROW_ORDEREDINPUT_H
