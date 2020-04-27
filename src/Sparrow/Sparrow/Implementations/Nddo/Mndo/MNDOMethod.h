/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MNDOMETHOD_H
#define SPARROW_MNDOMETHOD_H

#include <Utils/Scf/MethodInterfaces/ScfMethod.h>

namespace Scine {
namespace Utils {
enum class derivOrder;
}
namespace Sparrow {

namespace nddo {
class FockMatrix;
class NDDOInitializer;
class RawParametersContainer;
class OneElectronMatrix;
class TwoElectronMatrix;

class MNDOMethod : public Utils::ScfMethod {
 public:
  MNDOMethod();
  ~MNDOMethod() override;

  /*! Deprecated! Initialize the method from a structure and a parameter file. */
  void setStructure(const Utils::AtomCollection& atoms, std::string parameterPath);
  /*! Load the parameters from a file. */
  void readParameters(const std::string& parameterPath);
  /*! Save the parameters to a file. */
  void saveParameters(const std::string& fileName);

  NDDOInitializer& getInitializer() {
    return *mndoSettings_;
  }
  const NDDOInitializer& getInitializer() const {
    return *mndoSettings_;
  }

  /*! Get reference to the class for raw MNDO parameters. */
  RawParametersContainer& getRawParameters();
  /*! Get const reference to the class for raw MNDO parameters. */
  const RawParametersContainer& getRawParameters() const;

  const nddo::OneElectronMatrix& getOneElectronMatrix() const;
  const nddo::TwoElectronMatrix& getTwoElectronMatrix() const;

 private:
  std::shared_ptr<NDDOInitializer> mndoSettings_;
  std::shared_ptr<FockMatrix> mndoFock_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_MNDOMETHOD_H
