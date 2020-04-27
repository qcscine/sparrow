/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_PM6METHOD_H
#define SPARROW_PM6METHOD_H

#include <Utils/Scf/MethodInterfaces/ScfMethod.h>

namespace Scine {
namespace Sparrow {

namespace nddo {
class FockMatrix;
class NDDOInitializer;
class RawParametersContainer;
class OneElectronMatrix;
class TwoElectronMatrix;

class PM6Method : public Utils::ScfMethod {
 public:
  PM6Method();
  ~PM6Method() override;

  /*! Deprecated! Initialize the method from a structure and a parameter file. */
  void setStructure(const Utils::AtomCollection& atoms, std::string parameterPath);
  /*! Load the parameters from a file. */
  void readParameters(const std::string& parameterPath);
  /*! Save the parameters to a file. */
  void saveParameters(const std::string& fileName);

  const nddo::OneElectronMatrix& getOneElectronMatrix() const;
  const nddo::TwoElectronMatrix& getTwoElectronMatrix() const;

  NDDOInitializer& getInitializer() {
    return *pm6Settings_;
  }
  const NDDOInitializer& getInitializer() const {
    return *pm6Settings_;
  }

  /*! Get reference to the class for raw PM6 parameters. */
  RawParametersContainer& getRawParameters();
  /*! Get const reference to the class for raw PM6 parameters. */
  const RawParametersContainer& getRawParameters() const;

 private:
  std::shared_ptr<NDDOInitializer> pm6Settings_;
  std::shared_ptr<FockMatrix> pm6Fock_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_PM6METHOD_H
