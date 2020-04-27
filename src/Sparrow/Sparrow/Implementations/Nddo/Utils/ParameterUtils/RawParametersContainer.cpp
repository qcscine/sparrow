/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RawParametersContainer.h"
#include "RawParameters.h"
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Scf/MethodExceptions.h>
#include <algorithm>
#include <cassert>
#include <cereal/archives/xml.hpp>
#include <cereal/types/map.hpp>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Sparrow {

namespace nddo {

RawParametersContainer::RawParametersContainer(std::string path) : path_(path) {
  if (!path.empty()) {
    std::ifstream fs(path_);
    if (!fs)
      throw Utils::Methods::ParameterFileCannotBeOpenedException(path_);
    cereal::XMLInputArchive archive(fs);
    archive(atomics_, diatomics_);
  }
};

RawParametersContainer::~RawParametersContainer() = default;

RawParametersContainer::RawParametersContainer(RawParametersContainer&& rhs) noexcept = default;

RawParametersContainer& RawParametersContainer::operator=(RawParametersContainer&& rhs) = default;

void RawParametersContainer::writeParameterXMLFile(std::string fileName) {
  std::ofstream fs(fileName);
  cereal::XMLOutputArchive archive(fs);
  archive(atomics_, diatomics_);
}

bool RawParametersContainer::isAvailable(Utils::ElementType e) const {
  return atomics_.find(std::to_string(Utils::ElementInfo::Z(e))) != atomics_.end();
}

bool RawParametersContainer::isAvailable(Utils::ElementType e1, Utils::ElementType e2) const {
  std::string key(std::to_string(firstIndex(e1, e2)) + "-" + std::to_string(secondIndex(e1, e2)));
  return diatomics_.find(key) != diatomics_.end();
}

RawAtomicParameters& RawParametersContainer::getAtomicParameters(Utils::ElementType e) {
  assert(isAvailable(e) && "Parameters do not exist for e!");
  return atomics_.at(std::to_string(Utils::ElementInfo::Z(e)));
}

const RawAtomicParameters& RawParametersContainer::getAtomicParameters(Utils::ElementType e) const {
  assert(isAvailable(e) && "Parameters do not exist for e!");
  return atomics_.at(std::to_string(Utils::ElementInfo::Z(e)));
}

RawDiatomicParameters& RawParametersContainer::getDiatomicParameters(Utils::ElementType e1, Utils::ElementType e2) {
  assert(isAvailable(e1, e2) && "Parameters do not exist for e1-e2!");
  std::string key(std::to_string(Utils::ElementInfo::Z(e1)) + "-" + std::to_string(Utils::ElementInfo::Z(e2)));
  return diatomics_.at(key);
}

const RawDiatomicParameters& RawParametersContainer::getDiatomicParameters(Utils::ElementType e1, Utils::ElementType e2) const {
  assert(isAvailable(e1, e2) && "Parameters do not exist for e1-e2!");
  std::string key(std::to_string(Utils::ElementInfo::Z(e1)) + "-" + std::to_string(Utils::ElementInfo::Z(e2)));
  return diatomics_.at(key);
}

void RawParametersContainer::setAtomicParameters(Utils::ElementType e, const RawAtomicParameters& par) {
  atomics_[std::to_string(Utils::ElementInfo::Z(e))] = par;
}

void RawParametersContainer::setDiatomicParameters(Utils::ElementType e1, Utils::ElementType e2,
                                                   const RawDiatomicParameters& par) {
  std::string key(std::to_string(Utils::ElementInfo::Z(e1)) + "-" + std::to_string(Utils::ElementInfo::Z(e2)));
  diatomics_[key] = par;
}

unsigned int RawParametersContainer::index(Utils::ElementType e) const {
  return Utils::ElementInfo::Z(e);
}

unsigned int RawParametersContainer::firstIndex(Utils::ElementType e1, Utils::ElementType e2) const {
  auto i1 = Utils::ElementInfo::Z(e1);
  auto i2 = Utils::ElementInfo::Z(e2);

  return std::max(i1, i2);
}

unsigned int RawParametersContainer::secondIndex(Utils::ElementType e1, Utils::ElementType e2) const {
  auto i1 = Utils::ElementInfo::Z(e1);
  auto i2 = Utils::ElementInfo::Z(e2);

  return std::min(i1, i2);
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
