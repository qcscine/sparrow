/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Local2c2eIntegralCalculator.h"
#include "Sparrow/Implementations/Nddo/Utils/ParameterUtils/ChargeSeparationParameter.h"
#include "Sparrow/Implementations/Nddo/Utils/ParameterUtils/KlopmanParameter.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <cmath>

namespace Scine {
namespace Sparrow {

#ifdef REPRODUCE_PM6_ERRORS
#  define DONOTREPRODUCEERROR false // set to true in order not to reproduce mistakes of thiel1992a
#else
#  define DONOTREPRODUCEERROR true
#endif

using namespace Utils::AutomaticDifferentiation;

namespace nddo {
using namespace GeneralTypes;

namespace multipole {

template<Utils::DerivativeOrder O>
Value1DType<O> Local2c2eIntegralCalculator::getIntegral(twoElIntegral_t t1, twoElIntegral_t t2, double R,
                                                        const ChargeSeparationParameter& D1,
                                                        const ChargeSeparationParameter& d2,
                                                        const KlopmanParameter& rho1, const KlopmanParameter& rho2) {
  return getIntegral<O>(static_cast<int>(t1), static_cast<int>(t2), R, D1, d2, rho1, rho2);
}

template<Utils::DerivativeOrder O>
Value1DType<O> Local2c2eIntegralCalculator::getIntegral(int t1, int t2, double R, const ChargeSeparationParameter& D1,
                                                        const ChargeSeparationParameter& d2,
                                                        const KlopmanParameter& rho1, const KlopmanParameter& rho2) {
  static LocalTermArray terms = setUpTerms();
  auto sum = constant1D<O>(0);
  for (const auto& t : terms[t1][t2]) {
    double squaredRhos = (rho1.get(t.p1) + rho2.get(t.p2)) * (rho1.get(t.p1) + rho2.get(t.p2));
    double dist1 = t.p1 < MultipolePair::ss0 ? D1.get(t.p1) : 0;
    double dist2 = t.p2 < MultipolePair::ss0 ? d2.get(t.p2) : 0;
    sum += t.f * MultipoleMultipoleInteractionContainer::calculate<O>(t.m1, t.m2, R, dist1, dist2, squaredRhos);
  }
  return sum;
}

Local2c2eIntegralCalculator::LocalTermArray Local2c2eIntegralCalculator::setUpTerms() {
  LocalTermArray terms;

  for (int i = 0; i < 40; i++) {
    for (int j = 0; j < 40; j++) {
      // Charge configurations are ordered in such a way that the first 10 do not involve d basis functions
      bool hasD = !(i < 10 && j < 10);

      // List of pairs containing a prefactor and the charge configuration as a Multipole.
      auto l1 = getMultipoles(static_cast<twoElIntegral_t>(i), hasD);
      auto l2 = getMultipoles(static_cast<twoElIntegral_t>(j), hasD);

      auto op1 = separatePair(static_cast<twoElIntegral_t>(i));
      auto op2 = separatePair(static_cast<twoElIntegral_t>(j));
      // In integral of the form <o1a o1b |1/r| o2a o2b>
      orb_t o1a = op1.first;
      orb_t o1b = op1.second;
      orb_t o2a = op2.first;
      orb_t o2b = op2.second;
      // Orbital quantum numbers
      int l1a = orbitalQN(o1a);
      int l1b = orbitalQN(o1b);
      int l2a = orbitalQN(o2a);
      int l2b = orbitalQN(o2b);

      for (const auto& p1 : l1) {   // For single term of the multipolar expansion on first multipolar expansion
        for (const auto& p2 : l2) { // For single term of the multipolar expansion on second multipolar expansion
          // Get the charge configuration type of the specific term of the multipole expansion
          MultipolePair mp1 = pairType(l1a, l1b, LQuantumNumber(p1.second));
          MultipolePair mp2 = pairType(l2a, l2b, LQuantumNumber(p2.second));

          // [M_{20}, M_{20}] = [\~Q_{zx},\~Q_{zx}] + 3/4*[\~Q_{xy},\~Q_{xy}] is correct,
          // here [M_{20}, M_{20}] = [\~Q_{zx},\~Q_{zx}] - 1/4*[\~Q_{xy},\~Q_{xy}]
          // \~Q_{xy} = M22, \~Q_{zx} = Q_{zx} in this code snippet.
          if (hasD && p1.second == Multipole::M20 && p2.second == Multipole::M20 && !DONOTREPRODUCEERROR) { // Reproduce
                                                                                                            // error in
                                                                                                            // thiel1992a
            LocalTerm ta{}, tb{};
            ta.f = p1.first * p2.first;
            ta.m1 = Multipole::Qzx;
            ta.m2 = Multipole::Qzx;
            ta.p1 = mp1;
            ta.p2 = mp2;
            tb.f = p1.first * p2.first * (-0.25); // Instead of +0.75 as would be correct
            tb.m1 = Multipole::M22;
            tb.m2 = Multipole::M22;
            tb.p1 = mp1;
            tb.p2 = mp2;

            terms[i][j].push_back(ta);
            terms[i][j].push_back(tb);
          }
          else { // If the multipoles are not orthogonal, add them to the list of integrals to calculate.
            if (MQuantumNumber(p1.second) == MQuantumNumber(p2.second)) {
              LocalTerm t{};
              t.f = p1.first * p2.first;
              t.m1 = p1.second;
              t.m2 = p2.second;
              t.p1 = mp1;
              t.p2 = mp2;

              terms[i][j].push_back(t);
            }
          }
        }
      }
    }
  }

  return terms;
}

std::list<std::pair<double, Multipole>> Local2c2eIntegralCalculator::getMultipoles(twoElIntegral_t t, bool hasD) {
  std::list<std::pair<double, Multipole>> l;
  switch (t) {
    case twoElIntegral_t::s_s:
      l.emplace_back(1, Multipole::M00);
      break;
    case twoElIntegral_t::s_z:
      l.emplace_back(1, Multipole::M10);
      break;
    case twoElIntegral_t::s_x:
      l.emplace_back(1, Multipole::M11);
      break;
    case twoElIntegral_t::s_y:
      l.emplace_back(1, Multipole::M1m1);
      break;
    case twoElIntegral_t::x_z:
      l.emplace_back(1, Multipole::M21);
      break;
    case twoElIntegral_t::y_z:
      l.emplace_back(1, Multipole::M2m1);
      break;
    case twoElIntegral_t::x_y:
      l.emplace_back(1, Multipole::M2m2);
      break;

    case twoElIntegral_t::z_z:
      l.emplace_back(1, Multipole::M00);
      if (hasD)
        l.emplace_back(4.0 / 3.0, Multipole::M20);
      else
        l.emplace_back(1, Multipole::Qzz);
      break;

    case twoElIntegral_t::x_x:
      l.emplace_back(1, Multipole::M00);
      if (hasD) {
        l.emplace_back(-2.0 / 3.0, Multipole::M20);
        l.emplace_back(1, Multipole::M22);
      }
      else {
        l.emplace_back(1, Multipole::Qxx);
      }
      break;

    case twoElIntegral_t::y_y:
      l.emplace_back(1, Multipole::M00);
      if (hasD) {
        l.emplace_back(-2.0 / 3.0, Multipole::M20);
        l.emplace_back(-1, Multipole::M22);
      }
      else {
        l.emplace_back(1, Multipole::Qyy);
      }
      break;

    case twoElIntegral_t::s_z2:
      l.emplace_back(2 * (1 / std::sqrt(3.0)), Multipole::M20);
      break;
    case twoElIntegral_t::s_xz:
      l.emplace_back(1, Multipole::M21);
      break;
    case twoElIntegral_t::s_yz:
      l.emplace_back(1, Multipole::M2m1);
      break;
    case twoElIntegral_t::s_x2y2:
      l.emplace_back(1, Multipole::M22);
      break;
    case twoElIntegral_t::s_xy:
      l.emplace_back(1, Multipole::M2m2);
      break;

    case twoElIntegral_t::z_z2:
      l.emplace_back(2 * (1 / std::sqrt(3.0)), Multipole::M10);
      break;
    case twoElIntegral_t::z_xz:
      l.emplace_back(1, Multipole::M11);
      break;
    case twoElIntegral_t::z_yz:
      l.emplace_back(1, Multipole::M1m1);
      break;

    case twoElIntegral_t::x_z2:
      l.emplace_back(-(1 / std::sqrt(3.0)), Multipole::M11);
      break;
    case twoElIntegral_t::x_xz:
      l.emplace_back(1, Multipole::M10);
      break;
    case twoElIntegral_t::x_x2y2:
      l.emplace_back(1, Multipole::M11);
      break;
    case twoElIntegral_t::x_xy:
      l.emplace_back(1, Multipole::M1m1);
      break;

    case twoElIntegral_t::y_z2:
      l.emplace_back(-(1 / std::sqrt(3.0)), Multipole::M1m1);
      break;
    case twoElIntegral_t::y_yz:
      l.emplace_back(1, Multipole::M10);
      break;
    case twoElIntegral_t::y_x2y2:
      l.emplace_back(-1, Multipole::M1m1);
      break;
    case twoElIntegral_t::y_xy:
      l.emplace_back(1, Multipole::M11);
      break;

    case twoElIntegral_t::z2_z2:
      l.emplace_back(1, Multipole::M00);
      l.emplace_back(4.0 / 3.0, Multipole::M20);
      break;

    case twoElIntegral_t::z2_xz:
      l.emplace_back((1 / std::sqrt(3.0)), Multipole::M21);
      break;
    case twoElIntegral_t::z2_yz:
      l.emplace_back((1 / std::sqrt(3.0)), Multipole::M2m1);
      break;
    case twoElIntegral_t::z2_x2y2:
      l.emplace_back(-2 * (1 / std::sqrt(3.0)), Multipole::M22);
      break;
    case twoElIntegral_t::z2_xy:
      l.emplace_back(-2 * (1 / std::sqrt(3.0)), Multipole::M2m2);
      break;

    case twoElIntegral_t::xz_xz:
      l.emplace_back(1, Multipole::M00);
      l.emplace_back(2.0 / 3.0, Multipole::M20);
      l.emplace_back(1, Multipole::M22);
      break;

    case twoElIntegral_t::xz_yz:
      l.emplace_back(1, Multipole::M2m2);
      break;
    case twoElIntegral_t::xz_x2y2:
      l.emplace_back(1, Multipole::M21);
      break;
    case twoElIntegral_t::xz_xy:
      l.emplace_back(1, Multipole::M2m1);
      break;

    case twoElIntegral_t::yz_yz:
      l.emplace_back(1, Multipole::M00);
      l.emplace_back(2.0 / 3.0, Multipole::M20);
      l.emplace_back(-1, Multipole::M22);
      break;

    case twoElIntegral_t::yz_x2y2:
      l.emplace_back(-1, Multipole::M2m1);
      break;
    case twoElIntegral_t::yz_xy:
      l.emplace_back(1, Multipole::M21);
      break;

    case twoElIntegral_t::x2y2_x2y2:
      l.emplace_back(1, Multipole::M00);
      l.emplace_back(-4.0 / 3.0, Multipole::M20);
      break;

    case twoElIntegral_t::xy_xy:
      l.emplace_back(1, Multipole::M00);
      l.emplace_back(-4.0 / 3.0, Multipole::M20);
      break;
    default:
      break;
  }
  return l;
}

template Value1DType<Utils::DerivativeOrder::Zero> Local2c2eIntegralCalculator::getIntegral<Utils::DerivativeOrder::Zero>(
    twoElIntegral_t t1, twoElIntegral_t t2, double R, const ChargeSeparationParameter& D1,
    const ChargeSeparationParameter& d2, const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::DerivativeOrder::One> Local2c2eIntegralCalculator::getIntegral<Utils::DerivativeOrder::One>(
    twoElIntegral_t t1, twoElIntegral_t t2, double R, const ChargeSeparationParameter& D1,
    const ChargeSeparationParameter& d2, const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::DerivativeOrder::Two> Local2c2eIntegralCalculator::getIntegral<Utils::DerivativeOrder::Two>(
    twoElIntegral_t t1, twoElIntegral_t t2, double R, const ChargeSeparationParameter& D1,
    const ChargeSeparationParameter& d2, const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::DerivativeOrder::Zero> Local2c2eIntegralCalculator::getIntegral<Utils::DerivativeOrder::Zero>(
    int t1, int t2, double R, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& d2,
    const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::DerivativeOrder::One> Local2c2eIntegralCalculator::getIntegral<Utils::DerivativeOrder::One>(
    int t1, int t2, double R, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& d2,
    const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::DerivativeOrder::Two> Local2c2eIntegralCalculator::getIntegral<Utils::DerivativeOrder::Two>(
    int t1, int t2, double R, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& d2,
    const KlopmanParameter& rho1, const KlopmanParameter& rho2);

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
