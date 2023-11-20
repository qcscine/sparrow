/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Sparrow/Implementations/Nddo/Parameters.h"
#include "Utils/Geometry/ElementInfo.h"
#include "boost/filesystem.hpp"
#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#include <fstream>

namespace Scine {
namespace Sparrow {
namespace nddo {

/* Implementation of all archive templates */
template<typename Archive>
void Parameters::Atomic::Pack::Spd::serialize(Archive& archive) {
  archive(CEREAL_NVP(s), CEREAL_NVP(p), CEREAL_NVP(d));
}

template void Parameters::Atomic::Pack::Spd::serialize<cereal::JSONOutputArchive>(cereal::JSONOutputArchive&);
template void Parameters::Atomic::Pack::Spd::serialize<cereal::JSONInputArchive>(cereal::JSONInputArchive&);

template<typename Archive>
void Parameters::Atomic::Pack::serialize(Archive& archive) {
  archive(CEREAL_NVP(oneCenterEnergy), CEREAL_NVP(beta), CEREAL_NVP(orbitalExponent), CEREAL_NVP(internalExponent),
          CEREAL_NVP(gss), CEREAL_NVP(gpp), CEREAL_NVP(gsp), CEREAL_NVP(gp2), CEREAL_NVP(hsp), CEREAL_NVP(pcore),
          CEREAL_NVP(f0sd), CEREAL_NVP(g2sd), CEREAL_NVP(alpha));
}

template void Parameters::Atomic::Pack::serialize<cereal::JSONOutputArchive>(cereal::JSONOutputArchive&);
template void Parameters::Atomic::Pack::serialize<cereal::JSONInputArchive>(cereal::JSONInputArchive&);

template<typename Archive>
void Parameters::Atomic::GaussianRepulsion::serialize(Archive& archive) {
  archive(CEREAL_NVP(a), CEREAL_NVP(b), CEREAL_NVP(c));
}

template void Parameters::Atomic::GaussianRepulsion::serialize<cereal::JSONOutputArchive>(cereal::JSONOutputArchive&);
template void Parameters::Atomic::GaussianRepulsion::serialize<cereal::JSONInputArchive>(cereal::JSONInputArchive&);

template<typename Archive>
void Parameters::Atomic::serialize(Archive& archive) {
  archive(CEREAL_NVP(pack), CEREAL_NVP(gaussianRepulsion));
}

template void Parameters::Atomic::serialize<cereal::JSONOutputArchive>(cereal::JSONOutputArchive&);
template void Parameters::Atomic::serialize<cereal::JSONInputArchive>(cereal::JSONInputArchive&);

template<class Archive>
void Parameters::Diatomic::serialize(Archive& archive) {
  archive(CEREAL_NVP(exponent), CEREAL_NVP(factor));
}

template void Parameters::Diatomic::serialize<cereal::JSONOutputArchive>(cereal::JSONOutputArchive&);
template void Parameters::Diatomic::serialize<cereal::JSONInputArchive>(cereal::JSONInputArchive&);

Parameters::DiatomicKey Parameters::key(const int Z1, const int Z2) {
  if (Z1 > Z2) {
    return std::make_pair(Z1, Z2);
  }

  return std::make_pair(Z2, Z1);
}

Parameters::DiatomicKey Parameters::key(const Utils::ElementType a, const Utils::ElementType b) {
  return key(Utils::ElementInfo::Z(a), Utils::ElementInfo::Z(b));
}

Parameters Parameters::read(const std::string& filename) {
  if (!boost::filesystem::exists(filename)) {
    throw std::runtime_error("Parameter file to read does not exist");
  }

  std::ifstream infile(filename);
  cereal::JSONInputArchive archive(infile);
  Parameters result;
  archive(cereal::make_nvp("atomic", result.atomic), cereal::make_nvp("diatomic", result.diatomic));
  return result;
}

void Parameters::write(const std::string& filename) const {
  std::ofstream outfile(filename);
  cereal::JSONOutputArchive archive(outfile);
  archive(cereal::make_nvp("atomic", atomic), cereal::make_nvp("diatomic", diatomic));
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
