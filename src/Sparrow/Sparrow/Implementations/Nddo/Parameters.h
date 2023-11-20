/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SPARROW_NDDO_PARAMETERS_H
#define INCLUDE_SPARROW_NDDO_PARAMETERS_H

#include "Utils/Geometry/ElementTypes.h"
#include "boost/functional/hash.hpp"
#include <unordered_map>
#include <vector>

namespace Scine {
namespace Sparrow {
namespace nddo {

//! Nddo method parameters
struct Parameters {
  //!@name Member types {
  //!@{
  //! Atomic parameters
  struct Atomic {
    //!@name Member types
    //!@{
    /*! @brief Pack of implicitly constexpr part of atomic parameters
     *
     * @note In order to write conversion from constexpr parameters to runtime
     * parameters more easily, the trivial types whose basic operations are
     * implicitly constexpr are separated explicitly from those whose are not.
     */
    struct Pack {
      //!@name Member types
      //!@{
      struct Spd {
        double s = 0;
        double p = 0;
        double d = 0;

        template<class Archive>
        void serialize(Archive& archive);
      };
      //!@}

      //!@name Members
      //!@{
      //! One-center energy for each orbital in eV (uss, upp, udd)
      Spd oneCenterEnergy;
      //! Beta parameter for each orbital in eV (bs, bp, bd)
      Spd beta;
      //! Orbital exponent for each orbital in bohr^(-1) (zs, zp, zd)
      Spd orbitalExponent;
      //! Internal exponent for each orbital, needed for Slater-Condon in bohr^(-1) (zsn, zpn, zdn)
      Spd internalExponent;
      //! One-center Coulomb integral s-s (eV)
      double gss = 0;
      //! One-center Coulomb integral p-p (eV)
      double gpp = 0;
      //! One-center Coulomb integral s-p (eV)
      double gsp = 0;
      //! One-center Coulomb integral p-p' (eV)
      double gp2 = 0;
      //! One-center exchange integral s-s (eV)
      double hsp = 0;
      //! Klopman-Ohno term for core-core repulsion, sometimes different to pss0 (bohr)
      double pcore = 0;
      //! Slater-Condon parameter F0sd (eV)
      double f0sd = 0;
      //! Slater-Condon parameter G2sd (eV)
      double g2sd = 0;
      //! MNDO XXX
      double alpha = 0;
      //!@}

      template<class Archive>
      void serialize(Archive& archive);
    };

    struct GaussianRepulsion {
      double a = 0; /*!< Gaussian repulsion parameter a (-) */
      double b = 0; /*!< Gaussian repulsion parameter b (A^(-2)) */
      double c = 0; /*!< Gaussian repulsion parameter c (A) */

      template<class Archive>
      void serialize(Archive& archive);
    };
    //!@}

    //!@name Members
    //!@{
    //! Pack of constexpr parameters
    Pack pack;
    //! Gaussian repulsion parameters (usually 0-3)
    std::vector<GaussianRepulsion> gaussianRepulsion;
    //!@}

    template<class Archive>
    void serialize(Archive& archive);
  };

  //! Diatomic parameters
  struct Diatomic {
    //! Exponent of pairwise repulsion (A^(-1), for N-H / O-H / C-H A^(-2))
    double exponent = 0;
    //! Factor of pairwise repulsion
    double factor = 0;

    template<class Archive>
    void serialize(Archive& archive);
  };

  using DiatomicKey = std::pair<int, int>;
  //!@}

  //!@name Static functions
  //!@{
  //! Returns an ordered key to @p diatomic
  static DiatomicKey key(int Z1, int Z2);
  static DiatomicKey key(Utils::ElementType a, Utils::ElementType b);
  //! Read parameters from a JSON file
  static Parameters read(const std::string& filename);
  //!@}

  //!@name Observers
  //!@{
  //! Writes the parameters in JSON format to file
  void write(const std::string& filename) const;
  //!@}

  //!@name Members
  //!@{
  //! Map from atomic number to atomic parameters
  std::unordered_map<int, Atomic> atomic;
  //! Map from two atomic numbers to diatomic parameters
  std::unordered_map<DiatomicKey, Diatomic, boost::hash<DiatomicKey>> diatomic;
  //!@}
};

//! Default parameters for RM1
Parameters rm1();
//! Default parameters for AM1
Parameters am1();
//! Default parameters for PM3
Parameters pm3();
//! Default parameters for PM6
Parameters pm6();
//! Default parameters for MNDO
Parameters mndo();

} // namespace nddo
} // namespace Sparrow
} // namespace Scine

#endif
