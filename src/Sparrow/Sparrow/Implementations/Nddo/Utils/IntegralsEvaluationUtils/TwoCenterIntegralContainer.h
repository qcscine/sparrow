/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_TWOCENTERINTEGRALCONTAINER_H
#define SPARROW_TWOCENTERINTEGRALCONTAINER_H

#include "Global2c2eTerms.h"
#include <Utils/Typenames.h>
#include <memory>
#include <utility>
#include <vector>

namespace Scine {

namespace Sparrow {

namespace nddo {
class ElementParameters;
namespace multipole {
class Global2c2eMatrix;
}

/**
 * @brief This class contains smart pointers to two-center two-electron matrices for different atoms.
 *
 * This class stores for each atom pair a shared pointer to a Global2c2eMatrix, containing the ERIs between the
 * two centers.
 */
class TwoCenterIntegralContainer {
 public:
  using integralMatrix_t = std::shared_ptr<multipole::Global2c2eMatrix>;
  using Container = std::vector<std::vector<integralMatrix_t>>;

  /**
   * @brief constructor, give a reference to the positions, to the elements and to the element paramters, where the
   *        atomic orbital composition is.
   * @param elements vector specifying the elements of the molecule.
   * @param positions vector of the position in cartesian coordinates of the nuclei.
   * @param ep parameters of the elements. Contain, for istance, the AO composition of each element.
   */
  TwoCenterIntegralContainer(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                             const ElementParameters& ep);
  /**
   * @brief Initializes the Global2c2eMatrix for each atom pair. Their size depend on the element pair in question.
   */
  void initialize();
  /**
   * @brief Updated the Global2c2eMatrix for each atom pair.
   * @param order specify up to which derivative the integral has to be calculated.
   */
  void update(Utils::derivOrder order);

  /**
   * @brief sets a Global2c2eMatrix to correspond to a certain atom pair.
   * @param a index of the first atom.
   * @param b index of the second atom.
   * @param mat std::shared_ptr<multipole::Global2c2eMatrix> containing the ERIs corresponding to the atom pair.
   */
  void set(unsigned int a, unsigned int b, integralMatrix_t mat) {
    matrices_[a][b] = std::move(mat);
  }

  /**
   * @brief Getter for the ERIs corresponding to an atom pair.
   * @param a index of the first atom. Must be smaller than b.
   * @param b index of the second atom. Must be bigger than a.
   * @return a std::shared_ptr<Global2c2eMatrix> containing the ERIs corresponding to the atom pair.
   */
  integralMatrix_t get(unsigned int a, unsigned int b) const {
    // NB: TwoElectronMatrix designed such that always called with a<b.
    return matrices_[a][b];
  }

 private:
  // Initializes the matrix corresponding to a given atom pair.
  void initializePair(unsigned int i, unsigned int j);
  // Updates the matrix corresponding to a given atom pair.
  void updatePair(unsigned int i, unsigned int j, Utils::derivOrder order);

  multipole::Global2c2eTerms terms_;
  const ElementParameters& elementParameters_;
  Container matrices_;
  unsigned int nAtoms_;
  const Utils::ElementTypeCollection& elementTypes_;
  const Utils::PositionCollection& positions_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_TWOCENTERINTEGRALCONTAINER_H
