/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_RAWPM6PARAMETERS_H
#define SPARROW_RAWPM6PARAMETERS_H

#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace nddo {

/*!
 * Structure to store the atomic information present in the parameter file for PM6.
 */

struct gaussianRepulsionParameter {
  double grepa;
  double grepb;
  double grepc;

  template<class Archive>
  void serialize(Archive& archive) {
    archive(CEREAL_NVP(grepa), CEREAL_NVP(grepb), CEREAL_NVP(grepc));
  }
};

class RawAtomicParameters {
 public:
  double uss = 0,   /*!< One-center energy for s orbital (eV) */
      upp = 0,      /*!< One-center energy for s orbital (eV) */
      udd = 0;      /*!< One-center energy for s orbital (eV) */
  double bs = 0,    /*!< Beta parameter for s (eV) */
      bp = 0,       /*!< Beta parameter for p (eV) */
      bd = 0;       /*!< Beta parameter for d (eV) */
  double gss = 0,   /*!< One-center Coulomb integral s-s (eV) */
      gpp = 0,      /*!< One-center Coulomb integral p-p (eV) */
      gsp = 0,      /*!< One-center Coulomb integral s-p (eV) */
      gp2 = 0,      /*!< One-center Coulomb integral p-p' (eV) */
      hsp = 0;      /*!< One-center exchange integral s-s (eV) */
  double zs = 0,    /*!< Orbital exponent for s (bohr^(-1)) */
      zp = 0,       /*!< Orbital exponent for p (bohr^(-1)) */
      zd = 0;       /*!< Orbital exponent for d (bohr^(-1)) */
  double zsn = 0,   /*!< Internal exponent for s, needed for Slater-Condon parameters (bohr^(-1)) */
      zpn = 0,      /*!< Internal exponent for p, needed for Slater-Condon parameters (bohr^(-1)) */
      zdn = 0;      /*!< Internal exponent for d, needed for Slater-Condon parameters (bohr^(-1)) */
  double pcore = 0; /*!< Klopman-Ohno term for the core-core repulsion, sometimes unequal to pss0 (bohr) */
  double f0sd = 0,  /*!< Slater-Condon parameter F0sd (eV) */
      g2sd = 0;     /*!< Slater-Condon parameter G2sd (eV) */
  // double grepa = 0, /*!< Gaussian repulsion parameter a (-) */
  //    grepb = 0,    /*!< Gaussian repulsion parameter b (A^(-2)) */
  //    grepc = 0;    /*!< Gaussian repulsion parameter c (A) */
  double alpha = 0;                                                    /* MNDO XXX */
  std::vector<gaussianRepulsionParameter> gaussianRepulsionParameters; /*! Gaussian repulsion parameters */

  template<class Archive>
  void serialize(Archive& archive) {
    archive(CEREAL_NVP(uss), CEREAL_NVP(upp), CEREAL_NVP(udd), CEREAL_NVP(bs), CEREAL_NVP(bp), CEREAL_NVP(bd),
            CEREAL_NVP(gss), CEREAL_NVP(gpp), CEREAL_NVP(gsp), CEREAL_NVP(gp2), CEREAL_NVP(hsp), CEREAL_NVP(zs),
            CEREAL_NVP(zp), CEREAL_NVP(zd), CEREAL_NVP(zsn), CEREAL_NVP(zpn), CEREAL_NVP(zdn), CEREAL_NVP(pcore),
            CEREAL_NVP(f0sd), CEREAL_NVP(g2sd), CEREAL_NVP(alpha), CEREAL_NVP(gaussianRepulsionParameters));
  }
};

/*!
 * Structure to store the diatomic information present in the parameter file for PM6.
 */
class RawDiatomicParameters {
 public:
  double exponent = 0; /*!< Exponent of pairwise repulsion (A^(-1), for N-H / O-H / C-H A^(-2)) */
  double factor = 0;   /*!< Factor of pairwise repulsion (-) */

  template<class Archive>
  void serialize(Archive& archive) {
    archive(CEREAL_NVP(exponent), CEREAL_NVP(factor));
  }
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_RAWPM6PARAMETERS_H
