/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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

template<Utils::derivOrder O>
Value1DType<O> Local2c2eIntegralCalculator::getIntegral(twoElIntegral_t t1, twoElIntegral_t t2, double R,
                                                        const ChargeSeparationParameter& D1,
                                                        const ChargeSeparationParameter& d2,
                                                        const KlopmanParameter& rho1, const KlopmanParameter& rho2) {
  return getIntegral<O>(static_cast<int>(t1), static_cast<int>(t2), R, D1, d2, rho1, rho2);
}

template<Utils::derivOrder O>
Value1DType<O> Local2c2eIntegralCalculator::getIntegral(int t1, int t2, double R, const ChargeSeparationParameter& D1,
                                                        const ChargeSeparationParameter& d2,
                                                        const KlopmanParameter& rho1, const KlopmanParameter& rho2) {
  static LocalTermArray terms = setUpTerms();
  auto sum = constant1D<O>(0);
  for (const auto& t : terms[t1][t2]) {
    double squaredRhos = (rho1.get(t.p1) + rho2.get(t.p2)) * (rho1.get(t.p1) + rho2.get(t.p2));
    double dist1 = t.p1 < multipolePair_t::ss0 ? D1.get(t.p1) : 0;
    double dist2 = t.p2 < multipolePair_t::ss0 ? d2.get(t.p2) : 0;
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

      // List of pairs containing a prefactor and the charge configuration as a multipole_t.
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
          multipolePair_t mp1 = pairType(l1a, l1b, LQuantumNumber(p1.second));
          multipolePair_t mp2 = pairType(l2a, l2b, LQuantumNumber(p2.second));

          // [M_{20}, M_{20}] = [\~Q_{zx},\~Q_{zx}] + 3/4*[\~Q_{xy},\~Q_{xy}] is correct,
          // here [M_{20}, M_{20}] = [\~Q_{zx},\~Q_{zx}] - 1/4*[\~Q_{xy},\~Q_{xy}]
          // \~Q_{xy} = M22, \~Q_{zx} = Q_{zx} in this code snippet.
          if (hasD && p1.second == M20 && p2.second == M20 && !DONOTREPRODUCEERROR) { // Reproduce error in thiel1992a
            LocalTerm ta{}, tb{};
            ta.f = p1.first * p2.first;
            ta.m1 = Qzx;
            ta.m2 = Qzx;
            ta.p1 = mp1;
            ta.p2 = mp2;
            tb.f = p1.first * p2.first * (-0.25); // Instead of +0.75 as would be correct
            tb.m1 = M22;
            tb.m2 = M22;
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

std::list<std::pair<double, multipole_t>> Local2c2eIntegralCalculator::getMultipoles(twoElIntegral_t t, bool hasD) {
  std::list<std::pair<double, multipole_t>> l;
  switch (t) {
    case twoElIntegral_t::s_s:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      break;
    case twoElIntegral_t::s_z:
      l.push_back(std::make_pair<double, multipole_t>(1, M10));
      break;
    case twoElIntegral_t::s_x:
      l.push_back(std::make_pair<double, multipole_t>(1, M11));
      break;
    case twoElIntegral_t::s_y:
      l.push_back(std::make_pair<double, multipole_t>(1, M1m1));
      break;
    case twoElIntegral_t::x_z:
      l.push_back(std::make_pair<double, multipole_t>(1, M21));
      break;
    case twoElIntegral_t::y_z:
      l.push_back(std::make_pair<double, multipole_t>(1, M2m1));
      break;
    case twoElIntegral_t::x_y:
      l.push_back(std::make_pair<double, multipole_t>(1, M2m2));
      break;

    case twoElIntegral_t::z_z:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      if (hasD)
        l.push_back(std::make_pair<double, multipole_t>(4.0 / 3.0, M20));
      else
        l.push_back(std::make_pair<double, multipole_t>(1, Qzz));
      break;

    case twoElIntegral_t::x_x:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      if (hasD) {
        l.push_back(std::make_pair<double, multipole_t>(-2.0 / 3.0, M20));
        l.push_back(std::make_pair<double, multipole_t>(1, M22));
      }
      else {
        l.push_back(std::make_pair<double, multipole_t>(1, Qxx));
      }
      break;

    case twoElIntegral_t::y_y:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      if (hasD) {
        l.push_back(std::make_pair<double, multipole_t>(-2.0 / 3.0, M20));
        l.push_back(std::make_pair<double, multipole_t>(-1, M22));
      }
      else {
        l.push_back(std::make_pair<double, multipole_t>(1, Qyy));
      }
      break;

    case twoElIntegral_t::s_z2:
      l.push_back(std::make_pair<double, multipole_t>(2 * (1 / std::sqrt(3.0)), M20));
      break;
    case twoElIntegral_t::s_xz:
      l.push_back(std::make_pair<double, multipole_t>(1, M21));
      break;
    case twoElIntegral_t::s_yz:
      l.push_back(std::make_pair<double, multipole_t>(1, M2m1));
      break;
    case twoElIntegral_t::s_x2y2:
      l.push_back(std::make_pair<double, multipole_t>(1, M22));
      break;
    case twoElIntegral_t::s_xy:
      l.push_back(std::make_pair<double, multipole_t>(1, M2m2));
      break;

    case twoElIntegral_t::z_z2:
      l.push_back(std::make_pair<double, multipole_t>(2 * (1 / std::sqrt(3.0)), M10));
      break;
    case twoElIntegral_t::z_xz:
      l.push_back(std::make_pair<double, multipole_t>(1, M11));
      break;
    case twoElIntegral_t::z_yz:
      l.push_back(std::make_pair<double, multipole_t>(1, M1m1));
      break;

    case twoElIntegral_t::x_z2:
      l.push_back(std::make_pair<double, multipole_t>(-(1 / std::sqrt(3.0)), M11));
      break;
    case twoElIntegral_t::x_xz:
      l.push_back(std::make_pair<double, multipole_t>(1, M10));
      break;
    case twoElIntegral_t::x_x2y2:
      l.push_back(std::make_pair<double, multipole_t>(1, M11));
      break;
    case twoElIntegral_t::x_xy:
      l.push_back(std::make_pair<double, multipole_t>(1, M1m1));
      break;

    case twoElIntegral_t::y_z2:
      l.push_back(std::make_pair<double, multipole_t>(-(1 / std::sqrt(3.0)), M1m1));
      break;
    case twoElIntegral_t::y_yz:
      l.push_back(std::make_pair<double, multipole_t>(1, M10));
      break;
    case twoElIntegral_t::y_x2y2:
      l.push_back(std::make_pair<double, multipole_t>(-1, M1m1));
      break;
    case twoElIntegral_t::y_xy:
      l.push_back(std::make_pair<double, multipole_t>(1, M11));
      break;

    case twoElIntegral_t::z2_z2:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      l.push_back(std::make_pair<double, multipole_t>(4.0 / 3.0, M20));
      break;

    case twoElIntegral_t::z2_xz:
      l.push_back(std::make_pair<double, multipole_t>((1 / std::sqrt(3.0)), M21));
      break;
    case twoElIntegral_t::z2_yz:
      l.push_back(std::make_pair<double, multipole_t>((1 / std::sqrt(3.0)), M2m1));
      break;
    case twoElIntegral_t::z2_x2y2:
      l.push_back(std::make_pair<double, multipole_t>(-2 * (1 / std::sqrt(3.0)), M22));
      break;
    case twoElIntegral_t::z2_xy:
      l.push_back(std::make_pair<double, multipole_t>(-2 * (1 / std::sqrt(3.0)), M2m2));
      break;

    case twoElIntegral_t::xz_xz:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      l.push_back(std::make_pair<double, multipole_t>(2.0 / 3.0, M20));
      l.push_back(std::make_pair<double, multipole_t>(1, M22));
      break;

    case twoElIntegral_t::xz_yz:
      l.push_back(std::make_pair<double, multipole_t>(1, M2m2));
      break;
    case twoElIntegral_t::xz_x2y2:
      l.push_back(std::make_pair<double, multipole_t>(1, M21));
      break;
    case twoElIntegral_t::xz_xy:
      l.push_back(std::make_pair<double, multipole_t>(1, M2m1));
      break;

    case twoElIntegral_t::yz_yz:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      l.push_back(std::make_pair<double, multipole_t>(2.0 / 3.0, M20));
      l.push_back(std::make_pair<double, multipole_t>(-1, M22));
      break;

    case twoElIntegral_t::yz_x2y2:
      l.push_back(std::make_pair<double, multipole_t>(-1, M2m1));
      break;
    case twoElIntegral_t::yz_xy:
      l.push_back(std::make_pair<double, multipole_t>(1, M21));
      break;

    case twoElIntegral_t::x2y2_x2y2:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      l.push_back(std::make_pair<double, multipole_t>(-4.0 / 3.0, M20));
      break;

    case twoElIntegral_t::xy_xy:
      l.push_back(std::make_pair<double, multipole_t>(1, M00));
      l.push_back(std::make_pair<double, multipole_t>(-4.0 / 3.0, M20));
      break;
    default:
      break;
  }
  return l;
}

template Value1DType<Utils::derivOrder::zero> Local2c2eIntegralCalculator::getIntegral<Utils::derivOrder::zero>(
    twoElIntegral_t t1, twoElIntegral_t t2, double R, const ChargeSeparationParameter& D1,
    const ChargeSeparationParameter& d2, const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::derivOrder::one> Local2c2eIntegralCalculator::getIntegral<Utils::derivOrder::one>(
    twoElIntegral_t t1, twoElIntegral_t t2, double R, const ChargeSeparationParameter& D1,
    const ChargeSeparationParameter& d2, const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::derivOrder::two> Local2c2eIntegralCalculator::getIntegral<Utils::derivOrder::two>(
    twoElIntegral_t t1, twoElIntegral_t t2, double R, const ChargeSeparationParameter& D1,
    const ChargeSeparationParameter& d2, const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::derivOrder::zero> Local2c2eIntegralCalculator::getIntegral<Utils::derivOrder::zero>(
    int t1, int t2, double R, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& d2,
    const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::derivOrder::one> Local2c2eIntegralCalculator::getIntegral<Utils::derivOrder::one>(
    int t1, int t2, double R, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& d2,
    const KlopmanParameter& rho1, const KlopmanParameter& rho2);
template Value1DType<Utils::derivOrder::two> Local2c2eIntegralCalculator::getIntegral<Utils::derivOrder::two>(
    int t1, int t2, double R, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& d2,
    const KlopmanParameter& rho1, const KlopmanParameter& rho2);

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
