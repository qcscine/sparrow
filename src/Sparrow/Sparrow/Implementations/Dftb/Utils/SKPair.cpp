/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SKPair.h"
#include "SKAtom.h"
#include <Utils/IO/Regex.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/Technical/ScopedLocale.h>
#include <cmath>
#include <fstream>
#include <regex>

namespace Scine {
namespace Sparrow {

using std::exp;
using namespace Utils::AutomaticDifferentiation;

namespace dftb {

SKPair::SKPair(const std::string& atom1, const std::string& atom2, SKAtom* atomicParameters1, SKAtom* atomicParameters2,
               const std::string& path)
  : atomType1(atomicParameters1),
    atomType2(atomicParameters2),
    distFudge(1.0),
    extrC3(28),
    extrC4(28),
    extrC5(28),
    nInter(8),
    nInterRight(4),
    nInterLeft(4),
    deltaR(0.00001) {
  // Set integralIndexes in the order we need it (s orbitals first, then p, etc.)
  integralIndexes = {19, 9, 18, 8, 15, 5,  16, 6,  27, 23, 17, 7, 26, 22,
                     13, 3, 14, 4, 24, 20, 25, 21, 10, 0,  11, 1, 12, 2};

  // Read SKF file
  readSKFfile(atom1, atom2, atomicParameters1, atomicParameters2, path);

  // Set the number of integrals for the atom pair
  int nAOsA = atomicParameters1->getnAOs();
  int nAOsB = atomicParameters2->getnAOs();
  if (nAOsA + nAOsB == 2) // Only s orbitals
    nIntegrals = 2;
  else if (nAOsA + nAOsB == 5) // one atom only has s, the other one s and p
    nIntegrals = 4;
  else if (nAOsA + nAOsB == 8) // Both atoms have s and p orbitals
    nIntegrals = 10;
  else if (nAOsA + nAOsB == 10) // one has s, one has d
    nIntegrals = 14;
  else if (nAOsA + nAOsB == 13) // one has p, one has d
    nIntegrals = 22;
  else if (nAOsA + nAOsB == 18) // both have d
    nIntegrals = 28;
}

void SKPair::readSKFfile(const std::string& atom1, const std::string& atom2, SKAtom* atomicParameters1,
                         SKAtom* /*atomicParameters2*/, const std::string& path) {
  auto cLocale = Utils::ScopedLocale::cLocale();

  // construct filename and get it
  std::string filename = path + atom1 + "-" + atom2 + ".skf";
  std::ifstream fin(filename.c_str());

  // Output error if it doesn't exist
  if (!fin) {
    throw Utils::Methods::ParameterFileCannotBeOpenedException(filename);
  }

  // Variables used for going through the file
  std::string line;

  // Read first line
  getline(fin, line);
  if (line.find(".skf") != std::string::npos) { // Means that file points to another file
    fin.close();
    filename = path + line;
    fin.open(filename);
    getline(fin, line);
  }

  // The first line contains the grid distance and the number of grid points
  using namespace Utils::Regex;
  std::string firstLineRegexString = "^" + capturingFloatingPointNumber() + R"([,\s]+)" + capturingIntegerNumber();
  std::regex firstLineRegex(firstLineRegexString);
  std::smatch matches;
  std::regex_search(line, matches, firstLineRegex);
  gridDist = std::stod(matches[1]);
  nGridPoints = std::stoi(matches[2]);

  // If the atoms are identical, take supplementary information and put it into SKAtom
  if (atom1 == atom2) {
    getline(fin, line);
    std::string secondLineRegexString = "^" + capturingFloatingPointNumber();
    for (int i = 1; i < 10; ++i) {
      secondLineRegexString += R"([,\s]+)" + capturingFloatingPointNumber();
    }
    std::regex secondLineRegex(secondLineRegexString);
    std::regex_search(line, matches, secondLineRegex);
    atomicParameters1->setEnergies(std::stod(matches[3]), std::stod(matches[2]), std::stod(matches[1]));
    atomicParameters1->setHubbardParameter(std::stod(matches[7]), std::stod(matches[6]), std::stod(matches[5]));
    atomicParameters1->setOccupations(std::stoi(matches[10]), std::stoi(matches[9]), std::stoi(matches[8]));
  }

  // Read but ignore, so far, line containing mass and other parameters
  getline(fin, line);

  // Create array for the tabulated integrals
  M = std::vector<std::vector<double>>(28, std::vector<double>(nGridPoints + 50)); // Take 50 more to be sure there's
                                                                                   // enough space, in case the wrong
                                                                                   // number is given

  // Read the integrals
  nGridPoints--; // NB: -1 because some parameters set give nGridPoints value that is one too high...
  std::string oneEntryRegexString = R"((?:)" + capturingIntegerNumber() + R"(\*)?)" + capturingFloatingPointNumber();
  std::regex oneEntryRegex(oneEntryRegexString);

  for (int l = 0; l < nGridPoints; l++) {
    getline(fin, line);

    int p = 0;
    auto searchStart = line.cbegin();
    while (std::regex_search(searchStart, line.cend(), matches, oneEntryRegex)) {
      // first group is non-empty only if there is a multiplier, f.i. 5*0.0
      bool hasMultiplier = matches[1] != "";
      int multiplier = 1;
      double entry = std::stod(matches[2]);
      if (hasMultiplier) {
        multiplier = std::stoi(matches[1]);
      }

      for (int i = 0; i < multiplier; i++) {
        M[p++][l] = entry;
      }

      searchStart = matches[0].second;
    }
  }

  // Check first if it was really the end of integrals
  getline(fin, line);
  while (line[0] != 'S') {
    int p = 0;
    auto searchStart = line.cbegin();
    while (std::regex_search(searchStart, line.cend(), matches, oneEntryRegex)) {
      // first group is non-empty only if there is a multiplier, f.i. 5*0.0
      bool hasMultiplier = matches[1] != "";
      int multiplier = 1;
      double entry = std::stod(matches[2]);
      if (hasMultiplier) {
        multiplier = std::stoi(matches[1]);
      }

      for (int i = 0; i < multiplier; i++) {
        M[p++][nGridPoints] = entry;
      }

      searchStart = matches[0].second;
    }
    nGridPoints++;
    getline(fin, line);
  }

  rMax = gridDist * nGridPoints + distFudge;

  // Read spline parameters
  fin >> repulsion_.nSplineInts >> repulsion_.cutoff;
  fin >> repulsion_.a1 >> repulsion_.a2 >> repulsion_.a3;

  // Read spline (intervals)
  repulsion_.splineStart.resize(repulsion_.nSplineInts);
  repulsion_.splineEnd.resize(repulsion_.nSplineInts);
  repulsion_.c0.resize(repulsion_.nSplineInts);
  repulsion_.c1.resize(repulsion_.nSplineInts);
  repulsion_.c2.resize(repulsion_.nSplineInts);
  repulsion_.c3.resize(repulsion_.nSplineInts);
  for (int i = 0; i < repulsion_.nSplineInts; i++) {
    fin >> repulsion_.splineStart[i] >> repulsion_.splineEnd[i] >> repulsion_.c0[i] >> repulsion_.c1[i] >>
        repulsion_.c2[i] >> repulsion_.c3[i];
  }
  fin >> repulsion_.c4 >> repulsion_.c5;

  /*
  //Control:
  fin >> line;
  if(line != "<Documentation>")
    std::cout << "Error: not <Documentation>: " << line << std::endl;
  else
    std::cout << "Parameters read without error." << std::endl;
  */
}

void SKPair::precalculateGammaTerms() {
  // Precondition: not same atom type. In this case the gamma terms are not needed
  if (atomType1 == atomType2)
    return;

  double ta = atomType1->getHubbardParameter() * 3.2;
  double tb = atomType2->getHubbardParameter() * 3.2;
  double ta2 = ta * ta;
  double tb2 = tb * tb;

  g1a = (tb2 * tb2 * ta) / (2.0 * (ta2 - tb2) * (ta2 - tb2));
  g1b = (ta2 * ta2 * tb) / (2.0 * (tb2 - ta2) * (tb2 - ta2));
  g2a = tb2 * tb2 * (tb2 - 3 * ta2) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  g2b = ta2 * ta2 * (ta2 - 3 * tb2) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));

  dgab1a = (tb2 * tb2 * tb2 + 3.0 * ta2 * tb2 * tb2) / (2.0 * (ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  dgab2a = (12.0 * ta2 * ta * tb2 * tb2) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  dgab1b = (2.0 * ta2 * ta * tb2 * tb) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  dgab2b = (12.0 * ta2 * ta2 * tb2 * tb) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  dgba1a = (ta2 * ta2 * ta2 + 3.0 * tb2 * ta2 * ta2) / (2.0 * (tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));
  dgba2a = (12.0 * tb2 * tb * ta2 * ta2) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));
  dgba1b = (2.0 * tb2 * tb * ta2 * ta) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));
  dgba2b = (12.0 * tb2 * tb2 * ta2 * ta) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));

  dgadr = (tb2 * tb2 * tb2 - 3.0 * ta2 * tb2 * tb2) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  dgbdr = (ta2 * ta2 * ta2 - 3.0 * tb2 * ta2 * ta2) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));
}

void SKPair::getGammaTerms(double& t1a, double& t1b, double& t2a, double& t2b) {
  t1a = g1a;
  t1b = g1b;
  t2a = g2a;
  t2b = g2b;
}

void SKPair::getGammaDerTerms(double& dtab1a, double& dtab1b, double& dtab2a, double& dtab2b, double& dtba1a,
                              double& dtba1b, double& dtba2a, double& dtba2b, double& dtadr, double& dtbdr) {
  dtab1a = dgab1a;
  dtab1b = dgab1b;
  dtab2a = dgab2a;
  dtab2b = dgab2b;
  dtba1a = dgba1a;
  dtba1b = dgba1b;
  dtba2a = dgba2a;
  dtba2b = dgba2b;
  dtadr = dgadr;
  dtbdr = dgbdr;
}

void SKPair::complete(SKPair* p) {
  double factor = -1.0;
  if (p == this)
    factor = -1.0;
  for (int i = 0; i < nGridPoints; i++) {
    M[20][i] = factor * p->M[3][i];
    M[21][i] = factor * p->M[4][i];
    M[22][i] = p->M[7][i];
    M[23][i] = factor * p->M[8][i];
    M[24][i] = factor * p->M[13][i];
    M[25][i] = factor * p->M[14][i];
    M[26][i] = p->M[17][i];
    M[27][i] = factor * p->M[18][i];
  }
  // Precompute the coefficients for the extrapolation
  precompute5Extrapolation();
}

void SKPair::precompute5Extrapolation() {
  int indStart = nGridPoints - nInter;
  InterpolationValues<Utils::derivOrder::one> t1;
  InterpolationValues<Utils::derivOrder::one> t2;

  // Calculate interpolated values right after and right before the last given
  // point in parameter table; these values are used to compute the numerical
  // derivatives at the last point
  interpolate(t1, static_cast<double>(nGridPoints - 1) - deltaR, indStart); // Get values right before last given point
  interpolate(t2, static_cast<double>(nGridPoints - 1) + deltaR, indStart); // Get values right after last given point

  double yd0, yd1, yd2; // zero, first and second order derivatives
  double dx1, dx2;      // Temporary variables needed for the extrapolation

  for (int L = 0; L < nIntegrals; L++) {
    yd0 = M[integralIndexes[L]][nGridPoints - 1];
    yd1 = t1.derivIntegral[L].derivative();
    yd2 = (t2.derivIntegral[L].derivative() - t1.derivIntegral[L].derivative()) / (deltaR);

    // 5th order extrapolation (impose derivatives at last gridpoint, and zero value and derivatives at rMax)
    dx1 = -yd1 * (distFudge / gridDist);
    dx2 = yd2 * (distFudge / gridDist) * (distFudge / gridDist);
    extrC3[L] = 10 * yd0 - 4 * dx1 + 0.5 * dx2;
    extrC4[L] = -15 * yd0 + 7 * dx1 - dx2;
    extrC5[L] = 6 * yd0 - 3 * dx1 + 0.5 * dx2;
  }
}

template<Utils::derivOrder O>
Value1DType<O> SKPair::getRepulsion(double const& r) const {
  auto R = variableWithUnitDerivative<O>(r);
  if (r > repulsion_.cutoff)
    return constant1D<O>(0);
  if (r < repulsion_.splineStart[0])
    return exp(-repulsion_.a1 * R + repulsion_.a2) + repulsion_.a3;

  auto i = static_cast<int>((r - repulsion_.splineStart[0]) / (repulsion_.cutoff - repulsion_.splineStart[0]) *
                            repulsion_.nSplineInts);

  // If not the right bin, change it:
  if (repulsion_.splineStart[i] > r)
    while (repulsion_.splineStart[--i] > r) {
    }
  else if (repulsion_.splineEnd[i] < r)
    while (repulsion_.splineEnd[++i] < r) {
    }

  auto dr = R - repulsion_.splineStart[i];
  auto repulsion = repulsion_.c0[i] +
                   dr * (repulsion_.c1[i] +
                         dr * (repulsion_.c2[i] + dr * (repulsion_.c3[i] + (i == repulsion_.nSplineInts - 1
                                                                                ? dr * (repulsion_.c4 + dr * repulsion_.c5)
                                                                                : constant1D<O>(0.0)))));

  return repulsion;
}

template<Utils::derivOrder O>
int SKPair::getHS(double dist, InterpolationValues<O>& val) const {
  // If 'dist' too short or too large, return 0
  if (dist < gridDist || dist > rMax) {
    for (int L = 0; L < nIntegrals; L++)
      val.derivIntegral[L] = constant1D<O>(0.0);
    return 0;
  }

  // If 'dist' is in the extrapolation zone, take simple formula
  if (dist > nGridPoints * gridDist) {
    double dr = dist - rMax;
    auto xr = variableWithUnitDerivative<O>(-dr / distFudge);

    // 5th order extrapolation (impose derivatives at last gridpoint, and zero value and derivatives at rMax)
    for (int L = 0; L < nIntegrals; L++)
      val.derivIntegral[L] = ((extrC5[L] * xr + extrC4[L]) * xr + extrC3[L]) * xr * xr * xr;

    return 1;
  }

  // Else: do interpolation
  getHSIntegral(val, dist);

  return 1;
}

template<Utils::derivOrder O>
int SKPair::getHSIntegral(InterpolationValues<O>& val, double dist) const {
  double position = dist / gridDist - 1.0;
  int ind; // Current index during interpolation
  ind = static_cast<int>(position);

  int indStart; // first index for interpolation

  indStart = (ind > nGridPoints - nInterRight) ? nGridPoints - nInter : ((ind < nInterLeft) ? 0 : ind - nInterLeft + 1);
  interpolate(val, position, indStart);

  return 1;
}

template<Utils::derivOrder O>
void SKPair::interpolate(InterpolationValues<O>& val, double x, int start) const {
  /*
   * Taken from Numerical recipes
   * Adapted to treat simultaneously the different orbitals
   */
  double dift, dif = x - start;
  int ns = 0; // index of closest table entry
  for (int L = 0; L < nIntegrals; L++) {
    for (int i = 0; i < nInter; i++) {
      val.ya[L][i] = M[integralIndexes[L]][start + i];
      val.C[L][i] = constant1D<O>(M[integralIndexes[L]][start + i]);
      val.D[L][i] = constant1D<O>(M[integralIndexes[L]][start + i]);
    }
  }
  for (int i = 0; i < nInter; i++) {
    val.xa[i] = static_cast<double>(start + i);
    if ((dift = std::abs(x - val.xa[i])) < dif) {
      dif = dift;
      ns = i;
    }
  }
  for (int L = 0; L < nIntegrals; L++) {
    val.derivIntegral[L] = constant1D<O>(val.ya[L][ns]);
  }
  ns--;

  Value1DType<O> w, den;
  double ho, hp;
  for (int m = 1; m < nInter; m++) {
    for (int i = 0; i <= nInter - 1 - m; i++) {
      ho = val.xa[i] - x;
      hp = val.xa[i + m] - x;
      for (int L = 0; L < nIntegrals; L++) {
        w = val.C[L][i + 1] - val.D[L][i];
        den = w / (ho - hp);
        val.D[L][i] = getFromFull<O>(hp, -1.0 / gridDist, 0) * den;
        val.C[L][i] = getFromFull<O>(ho, -1.0 / gridDist, 0) * den;
      }
    }

    for (int L = 0; L < nIntegrals; L++) {
      val.derivIntegral[L] += (2 * (ns + 1) < (nInter - m)) ? val.C[L][ns + 1] : val.D[L][ns];
    }
    if (!(2 * (ns + 1) < (nInter - m)))
      ns--;
  }
}

template int SKPair::getHS<Utils::derivOrder::zero>(double, InterpolationValues<Utils::derivOrder::zero>&) const;
template int SKPair::getHS<Utils::derivOrder::one>(double, InterpolationValues<Utils::derivOrder::one>&) const;
template int SKPair::getHS<Utils::derivOrder::two>(double, InterpolationValues<Utils::derivOrder::two>&) const;
template int SKPair::getHSIntegral<Utils::derivOrder::zero>(InterpolationValues<Utils::derivOrder::zero>&, double) const;
template int SKPair::getHSIntegral<Utils::derivOrder::one>(InterpolationValues<Utils::derivOrder::one>&, double) const;
template int SKPair::getHSIntegral<Utils::derivOrder::two>(InterpolationValues<Utils::derivOrder::two>&, double) const;
template void SKPair::interpolate<Utils::derivOrder::zero>(InterpolationValues<Utils::derivOrder::zero>&, double, int) const;
template void SKPair::interpolate<Utils::derivOrder::one>(InterpolationValues<Utils::derivOrder::one>&, double, int) const;
template void SKPair::interpolate<Utils::derivOrder::two>(InterpolationValues<Utils::derivOrder::two>&, double, int) const;
template double SKPair::getRepulsion<Utils::derivOrder::zero>(double const& r) const;
template First1D SKPair::getRepulsion<Utils::derivOrder::one>(double const& r) const;
template Second1D SKPair::getRepulsion<Utils::derivOrder::two>(double const& r) const;

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
