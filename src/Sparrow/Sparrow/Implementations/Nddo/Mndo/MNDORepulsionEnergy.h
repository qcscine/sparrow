/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MNDOREPULSIONENERGY_H
#define SPARROW_MNDOREPULSIONENERGY_H

#include "MNDOPairwiseRepulsion.h"
#include <Utils/Scf/MethodInterfaces/RepulsionCalculator.h>
#include <memory>
#include <vector>

namespace Scine {
namespace Utils {
enum class derivOrder;
}
namespace Sparrow {

namespace nddo {
class ElementParameters;

/**
 * @brief This class sums up the core-core repulsion energies and the corresponding derivatives with respect to
 *        the nuclear cartesian coordinate between all pairs of cores.
 * It inherits from Utils::RepulsionCalculator in order for it to work with the LCAO/ScfMethod polymorphic system.
 */
class MNDORepulsionEnergy : public Utils::RepulsionCalculator {
 public:
  using pairRepulsion_t = std::unique_ptr<MNDOPairwiseRepulsion>;
  using Container = std::vector<std::vector<pairRepulsion_t>>;

  //! @brief Constructor.
  MNDORepulsionEnergy(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                      const ElementParameters& elementParameters);
  //! @brief Overrides virtual base class desctructor with default implementation.
  ~MNDORepulsionEnergy() override;

  //! @brief Initializes the core-core repulsion pairs
  void initialize() override;
  //! @brief Starts the calculation of the core-core repulsion up to the \param order derivative order.
  void calculateRepulsion(Utils::derivOrder order) override;
  //! @brief Sums up all the single core-core contributions to return the overall core-core repulsion energy.
  double getRepulsionEnergy() const override;
  //! Functions calculating the core-core derivative contributions up to the corresponding derivative order.
  void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const override;
  void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const override;
  void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const override;

 private:
  template<Utils::derivativeType O>
  void addRepulsionDerivativesImpl(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const;

  void calculatePairRepulsion(int i, int j, Utils::derivOrder order);
  void initializePair(int i, int j);

  const ElementParameters& elementParameters_;
  Container rep_;
  int nAtoms_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_REPULSIONENERGY_H
