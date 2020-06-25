/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_RAWPM6PARAMETERSCONTAINER
#define SPARROW_RAWPM6PARAMETERSCONTAINER

#include <Utils/Geometry/ElementTypes.h>
#include <array>
#include <map>
#include <memory>
#include <string>

namespace Scine {
namespace Sparrow {

namespace nddo {
class RawAtomicParameters;
class RawDiatomicParameters;

/*!
 * Class containing the PM6 parameters in a raw form.
 */
class RawParametersContainer {
 public:
  RawParametersContainer(std::string path = {});
  ~RawParametersContainer();
  RawParametersContainer(RawParametersContainer&& rhs) noexcept;
  RawParametersContainer& operator=(RawParametersContainer&& rhs);

  /*! Verify whether atomic parameters for e are available. */
  bool isAvailable(Utils::ElementType e) const;
  /*! Verify whether diatomic parameters for the pair e1 - e2 are available. */
  bool isAvailable(Utils::ElementType e1, Utils::ElementType e2) const;

  /*! Returns reference to already existing raw atomic parameters. */
  RawAtomicParameters& getAtomicParameters(Utils::ElementType e);
  /*! Returns const reference to already existing raw atomic parameters. */
  const RawAtomicParameters& getAtomicParameters(Utils::ElementType e) const;
  /*! Returns reference to already existing raw diatomic parameters. */
  RawDiatomicParameters& getDiatomicParameters(Utils::ElementType e1, Utils::ElementType e2);
  /*! Returns const reference to already existing raw diatomic parameters. */
  const RawDiatomicParameters& getDiatomicParameters(Utils::ElementType e1, Utils::ElementType e2) const;

  /*! Set atomic parameters for element type e. */
  void setAtomicParameters(Utils::ElementType e, const RawAtomicParameters& par);
  /*! Set diatomic parameters for element type pair e1 - e2. */
  void setDiatomicParameters(Utils::ElementType e1, Utils::ElementType e2, const RawDiatomicParameters& par);
  /**
   * @brief Writes the current parameters into an XML file.
   * @param fileName The name and path to the file to be written.
   */
  void writeParameterXMLFile(std::string fileName);

 private:
  unsigned int index(Utils::ElementType e) const;
  unsigned int firstIndex(Utils::ElementType e1, Utils::ElementType e2) const;
  unsigned int secondIndex(Utils::ElementType e1, Utils::ElementType e2) const;

  std::string path_;
  std::map<std::string, RawAtomicParameters> atomics_;
  std::map<std::string, RawDiatomicParameters> diatomics_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // RAWPM6PARAMETERSCONTAINER
