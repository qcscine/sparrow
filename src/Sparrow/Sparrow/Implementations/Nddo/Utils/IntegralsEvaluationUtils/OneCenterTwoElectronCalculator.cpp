/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OneCenterTwoElectronCalculator.h"
#include "GeneralTypes.h"
#include <cmath>

namespace Scine {
namespace Sparrow {

namespace nddo {
using namespace GeneralTypes;

bool OneCenterTwoElectronCalculator::indexesSet_ = false;
int OneCenterTwoElectronCalculator::index[9][9][9][9];
OneCenterTwoElectronIntegralExpression OneCenterTwoElectronCalculator::expr[58];

OneCenterTwoElectronCalculator::OneCenterTwoElectronCalculator() {
  if (!indexesSet_)
    setIndexes();
}

double OneCenterTwoElectronCalculator::getIntegral(int ind, const SlaterCondonParameters* p) {
  return expr[ind].result(p);
}

void OneCenterTwoElectronCalculator::setIndexes() {
  indexesSet_ = true;
  // Initialize all values of index to 100.
  for (auto& index1 : index) {
    for (auto& index2 : index1) {
      for (auto& index3 : index2) {
        for (auto& index4 : index3) {
          index4 = 100;
        }
      }
    }
  }

  setUniqueIndexes();

  // NB: STRANGE ERROR with Visual Studio Compiler
  // If I use 'i' directly instead of 'a', there is a compiler error (unknown)
  int i = 0;
  for (int a = 0; a < 9; a++) {
    i = a;
    for (int j = 0; j < 9; j++) {
      for (int k = 0; k < 9; k++) {
        for (int l = 0; l < 9; l++) {
          if (index[i][j][k][l] != 100) {
            setIndex(i, j, k, l, &index[j][i][k][l]);
            setIndex(i, j, k, l, &index[i][j][l][k]);
            setIndex(i, j, k, l, &index[j][i][l][k]);
            setIndex(i, j, k, l, &index[k][l][i][j]);
            setIndex(i, j, k, l, &index[k][l][j][i]);
            setIndex(i, j, k, l, &index[l][k][i][j]);
            setIndex(i, j, k, l, &index[l][k][j][i]);
          }
        }
      }
    }
  }
}

void OneCenterTwoElectronCalculator::setIndex(int i, int j, int k, int l, int* p) {
  *p = getIndex(i, j, k, l);
}

int OneCenterTwoElectronCalculator::getIndex(int i, int j, int k, int l) {
  return index[i][j][k][l];
}

void OneCenterTwoElectronCalculator::setUniqueIndexes() {
  index[orb_t::s][orb_t::s][orb_t::s][orb_t::s] = 0;
  expr[0] = OneCenterTwoElectronIntegralExpression(1, F0ss);

  index[orb_t::s][orb_t::s][orb_t::z][orb_t::z] = 1;
  index[orb_t::s][orb_t::s][orb_t::x][orb_t::x] = 1;
  index[orb_t::s][orb_t::s][orb_t::y][orb_t::y] = 1;
  expr[1] = OneCenterTwoElectronIntegralExpression(1, F0sp);

  index[orb_t::s][orb_t::z][orb_t::s][orb_t::z] = 2;
  index[orb_t::s][orb_t::x][orb_t::s][orb_t::x] = 2;
  index[orb_t::s][orb_t::y][orb_t::s][orb_t::y] = 2;
  expr[2] = OneCenterTwoElectronIntegralExpression(1.0 / 3.0, G1sp);

  index[orb_t::z][orb_t::z][orb_t::z][orb_t::z] = 3;
  index[orb_t::x][orb_t::x][orb_t::x][orb_t::x] = 3;
  index[orb_t::y][orb_t::y][orb_t::y][orb_t::y] = 3;
  expr[3] = OneCenterTwoElectronIntegralExpression(1, F0pp, 4.0 / 25.0, F2pp);

  index[orb_t::x][orb_t::x][orb_t::z][orb_t::z] = 4;
  index[orb_t::y][orb_t::y][orb_t::z][orb_t::z] = 4;
  index[orb_t::x][orb_t::x][orb_t::y][orb_t::y] = 4;
  expr[4] = OneCenterTwoElectronIntegralExpression(1, F0pp, -2.0 / 25.0, F2pp);

  index[orb_t::x][orb_t::z][orb_t::x][orb_t::z] = 5;
  index[orb_t::y][orb_t::z][orb_t::y][orb_t::z] = 5;
  index[orb_t::x][orb_t::y][orb_t::x][orb_t::y] = 5;
  expr[5] = OneCenterTwoElectronIntegralExpression(3.0 / 25.0, F2pp);

  index[orb_t::x][orb_t::z2][orb_t::x][orb_t::z2] = 6;
  index[orb_t::y][orb_t::z2][orb_t::y][orb_t::z2] = 6;
  expr[6] = OneCenterTwoElectronIntegralExpression(1.0 / 15.0, G1pd, 18.0 / 245.0, G3pd);

  index[orb_t::x][orb_t::z2][orb_t::x][orb_t::x2y2] = 7;
  index[orb_t::x][orb_t::z2][orb_t::y][orb_t::xy] = 7;
  index[orb_t::x][orb_t::xy][orb_t::y][orb_t::z2] = 7;
  expr[7] = OneCenterTwoElectronIntegralExpression(-sqrt(3) / 15.0, G1pd, -sqrt(27.0) / 245.0, G3pd);

  index[orb_t::x][orb_t::x2y2][orb_t::y][orb_t::xy] = 8;
  expr[8] = OneCenterTwoElectronIntegralExpression(1.0 / 5.0, G1pd, -21.0 / 245.0, G3pd);

  index[orb_t::x][orb_t::xy][orb_t::y][orb_t::x2y2] = 9;
  expr[9] = OneCenterTwoElectronIntegralExpression(-1.0 / 5.0, G1pd, 21.0 / 245.0, G3pd);

  index[orb_t::x][orb_t::x][orb_t::z2][orb_t::z2] = 10;
  index[orb_t::y][orb_t::y][orb_t::z2][orb_t::z2] = 10;
  expr[10] = OneCenterTwoElectronIntegralExpression(1, F0pd, -2.0 / 35.0, F2pd);

  index[orb_t::x][orb_t::x][orb_t::z2][orb_t::x2y2] = 11;
  index[orb_t::x][orb_t::y][orb_t::z2][orb_t::xy] = 11;
  expr[11] = OneCenterTwoElectronIntegralExpression(-sqrt(12.0) / 35.0, F2pd);

  index[orb_t::y][orb_t::z2][orb_t::y][orb_t::x2y2] = 12;
  expr[12] = OneCenterTwoElectronIntegralExpression(sqrt(3.0) / 15.0, G1pd, sqrt(27) / 245.0, G3pd);

  index[orb_t::z][orb_t::z2][orb_t::z][orb_t::z2] = 13;
  expr[13] = OneCenterTwoElectronIntegralExpression(4.0 / 15.0, G1pd, 27.0 / 245.0, G3pd);

  index[orb_t::z][orb_t::x2y2][orb_t::z][orb_t::x2y2] = 14;
  index[orb_t::z][orb_t::xy][orb_t::z][orb_t::xy] = 14;
  index[orb_t::x][orb_t::yz][orb_t::x][orb_t::yz] = 14;
  index[orb_t::y][orb_t::xz][orb_t::y][orb_t::xz] = 14;
  index[orb_t::x][orb_t::xz][orb_t::z][orb_t::x2y2] = 14;
  index[orb_t::x][orb_t::yz][orb_t::z][orb_t::xy] = 14;
  index[orb_t::y][orb_t::xz][orb_t::z][orb_t::xy] = 14;
  index[orb_t::x][orb_t::yz][orb_t::y][orb_t::xz] = 14;
  expr[14] = OneCenterTwoElectronIntegralExpression(3.0 / 49.0, G3pd);

  index[orb_t::z][orb_t::xz][orb_t::z][orb_t::xz] = 15;
  index[orb_t::z][orb_t::yz][orb_t::z][orb_t::yz] = 15;
  index[orb_t::x][orb_t::xy][orb_t::x][orb_t::xy] = 15;
  index[orb_t::x][orb_t::x2y2][orb_t::x][orb_t::x2y2] = 15;
  index[orb_t::x][orb_t::xz][orb_t::x][orb_t::xz] = 15;
  index[orb_t::y][orb_t::xy][orb_t::y][orb_t::xy] = 15;
  index[orb_t::y][orb_t::yz][orb_t::y][orb_t::yz] = 15;
  index[orb_t::y][orb_t::x2y2][orb_t::y][orb_t::x2y2] = 15;
  expr[15] = OneCenterTwoElectronIntegralExpression(1.0 / 5.0, G1pd, 24.0 / 245.0, G3pd);

  index[orb_t::x][orb_t::xz][orb_t::z][orb_t::z2] = 16;
  index[orb_t::y][orb_t::yz][orb_t::z][orb_t::z2] = 16;
  expr[16] = OneCenterTwoElectronIntegralExpression(sqrt(12.0) / 15.0, G1pd, -sqrt(243) / 245.0, G3pd);

  index[orb_t::x][orb_t::x2y2][orb_t::z][orb_t::xz] = 17;
  index[orb_t::x][orb_t::xy][orb_t::z][orb_t::yz] = 17;
  index[orb_t::y][orb_t::xy][orb_t::z][orb_t::xz] = 17;
  index[orb_t::x][orb_t::xz][orb_t::y][orb_t::yz] = 17;
  expr[17] = OneCenterTwoElectronIntegralExpression(1.0 / 5.0, G1pd, -6.0 / 245.0, G3pd);

  index[orb_t::y][orb_t::yz][orb_t::z][orb_t::x2y2] = 18;
  expr[18] = OneCenterTwoElectronIntegralExpression(-3.0 / 49.0, G3pd);

  index[orb_t::y][orb_t::z2][orb_t::z][orb_t::yz] = 19;
  index[orb_t::x][orb_t::z2][orb_t::z][orb_t::xz] = 19;
  expr[19] = OneCenterTwoElectronIntegralExpression(-sqrt(3) / 15.0, G1pd, sqrt(432) / 245.0, G3pd);

  index[orb_t::y][orb_t::x2y2][orb_t::z][orb_t::yz] = 20;
  expr[20] = OneCenterTwoElectronIntegralExpression(-1.0 / 5.0, G1pd, 6.0 / 245.0, G3pd);

  index[orb_t::z][orb_t::z][orb_t::z2][orb_t::z2] = 21;
  expr[21] = OneCenterTwoElectronIntegralExpression(1, F0pd, 4.0 / 35.0, F2pd);

  index[orb_t::x][orb_t::z][orb_t::z2][orb_t::xz] = 22;
  index[orb_t::y][orb_t::z][orb_t::z2][orb_t::yz] = 22;
  expr[22] = OneCenterTwoElectronIntegralExpression(sqrt(3) / 35.0, F2pd);

  index[orb_t::z][orb_t::z][orb_t::x2y2][orb_t::x2y2] = 23;
  index[orb_t::z][orb_t::z][orb_t::xy][orb_t::xy] = 23;
  index[orb_t::x][orb_t::x][orb_t::yz][orb_t::yz] = 23;
  index[orb_t::y][orb_t::y][orb_t::xz][orb_t::xz] = 23;
  expr[23] = OneCenterTwoElectronIntegralExpression(1, F0pd, -4.0 / 35.0, F2pd);

  index[orb_t::x][orb_t::z][orb_t::x2y2][orb_t::xz] = 24;
  index[orb_t::x][orb_t::z][orb_t::xy][orb_t::yz] = 24;
  index[orb_t::y][orb_t::z][orb_t::xy][orb_t::xz] = 24;
  index[orb_t::x][orb_t::y][orb_t::xz][orb_t::yz] = 24;
  expr[24] = OneCenterTwoElectronIntegralExpression(3.0 / 35.0, F2pd);

  index[orb_t::y][orb_t::z][orb_t::x2y2][orb_t::yz] = 25;
  expr[25] = OneCenterTwoElectronIntegralExpression(-3.0 / 35.0, F2pd);

  index[orb_t::z][orb_t::z][orb_t::xz][orb_t::xz] = 26;
  index[orb_t::z][orb_t::z][orb_t::yz][orb_t::yz] = 26;
  index[orb_t::x][orb_t::x][orb_t::xy][orb_t::xy] = 26;
  index[orb_t::x][orb_t::x][orb_t::x2y2][orb_t::x2y2] = 26;
  index[orb_t::x][orb_t::x][orb_t::xz][orb_t::xz] = 26;
  index[orb_t::y][orb_t::y][orb_t::xy][orb_t::xy] = 26;
  index[orb_t::y][orb_t::y][orb_t::yz][orb_t::yz] = 26;
  index[orb_t::y][orb_t::y][orb_t::x2y2][orb_t::x2y2] = 26;
  expr[26] = OneCenterTwoElectronIntegralExpression(1, F0pd, 2.0 / 35.0, F2pd);

  index[orb_t::y][orb_t::y][orb_t::z2][orb_t::x2y2] = 27;
  expr[27] = OneCenterTwoElectronIntegralExpression(sqrt(12) / 35.0, F2pd);

  index[orb_t::s][orb_t::z2][orb_t::s][orb_t::z2] = 28;
  index[orb_t::s][orb_t::x2y2][orb_t::s][orb_t::x2y2] = 28;
  index[orb_t::s][orb_t::xy][orb_t::s][orb_t::xy] = 28;
  index[orb_t::s][orb_t::xz][orb_t::s][orb_t::xz] = 28;
  index[orb_t::s][orb_t::yz][orb_t::s][orb_t::yz] = 28;
  expr[28] = OneCenterTwoElectronIntegralExpression(1.0 / 5.0, G2sd);

  index[orb_t::s][orb_t::s][orb_t::z2][orb_t::z2] = 29;
  index[orb_t::s][orb_t::s][orb_t::x2y2][orb_t::x2y2] = 29;
  index[orb_t::s][orb_t::s][orb_t::xy][orb_t::xy] = 29;
  index[orb_t::s][orb_t::s][orb_t::xz][orb_t::xz] = 29;
  index[orb_t::s][orb_t::s][orb_t::yz][orb_t::yz] = 29;
  expr[29] = OneCenterTwoElectronIntegralExpression(1, F0sd);

  index[orb_t::z2][orb_t::z2][orb_t::z2][orb_t::z2] = 30;
  index[orb_t::x2y2][orb_t::x2y2][orb_t::x2y2][orb_t::x2y2] = 30;
  index[orb_t::xy][orb_t::xy][orb_t::xy][orb_t::xy] = 30;
  index[orb_t::xz][orb_t::xz][orb_t::xz][orb_t::xz] = 30;
  index[orb_t::yz][orb_t::yz][orb_t::yz][orb_t::yz] = 30;
  expr[30] = OneCenterTwoElectronIntegralExpression(1, F0dd, 4.0 / 49.0, F2dd, 36.0 / 441.0, F4dd);

  index[orb_t::z2][orb_t::x2y2][orb_t::z2][orb_t::x2y2] = 31;
  index[orb_t::z2][orb_t::xy][orb_t::z2][orb_t::xy] = 31;
  expr[31] = OneCenterTwoElectronIntegralExpression(4.0 / 49.0, F2dd, 15.0 / 441.0, F4dd);

  index[orb_t::z2][orb_t::xz][orb_t::z2][orb_t::xz] = 32;
  index[orb_t::z2][orb_t::yz][orb_t::z2][orb_t::yz] = 32;
  expr[32] = OneCenterTwoElectronIntegralExpression(1.0 / 49.0, F2dd, 30.0 / 441.0, F4dd);

  index[orb_t::z2][orb_t::z2][orb_t::x2y2][orb_t::x2y2] = 33;
  index[orb_t::z2][orb_t::z2][orb_t::xy][orb_t::xy] = 33;
  expr[33] = OneCenterTwoElectronIntegralExpression(1, F0dd, -4.0 / 49.0, F2dd, 6.0 / 441.0, F4dd);

  index[orb_t::z2][orb_t::xz][orb_t::x2y2][orb_t::xz] = 34;
  index[orb_t::z2][orb_t::xz][orb_t::xy][orb_t::yz] = 34;
  index[orb_t::z2][orb_t::yz][orb_t::xy][orb_t::xz] = 34;
  expr[34] = OneCenterTwoElectronIntegralExpression(sqrt(3) / 49.0, F2dd, -sqrt(75) / 441.0, F4dd);

  index[orb_t::z2][orb_t::yz][orb_t::x2y2][orb_t::yz] = 35;
  expr[35] = OneCenterTwoElectronIntegralExpression(-sqrt(3) / 49.0, F2dd, sqrt(75) / 441.0, F4dd);

  index[orb_t::z2][orb_t::z2][orb_t::xz][orb_t::xz] = 36;
  index[orb_t::z2][orb_t::z2][orb_t::yz][orb_t::yz] = 36;
  expr[36] = OneCenterTwoElectronIntegralExpression(1, F0dd, 2.0 / 49.0, F2dd, -24.0 / 441.0, F4dd);

  index[orb_t::z2][orb_t::x2y2][orb_t::xz][orb_t::xz] = 37;
  index[orb_t::z2][orb_t::xy][orb_t::xz][orb_t::yz] = 37;
  expr[37] = OneCenterTwoElectronIntegralExpression(-sqrt(12) / 49.0, F2dd, sqrt(300) / 441.0, F4dd);

  index[orb_t::z2][orb_t::x2y2][orb_t::yz][orb_t::yz] = 38;
  expr[38] = OneCenterTwoElectronIntegralExpression(sqrt(12) / 49.0, F2dd, -sqrt(300) / 441.0, F4dd);

  index[orb_t::x2y2][orb_t::xy][orb_t::x2y2][orb_t::xy] = 39;
  expr[39] = OneCenterTwoElectronIntegralExpression(35.0 / 441.0, F4dd);

  index[orb_t::xy][orb_t::xz][orb_t::xy][orb_t::xz] = 40;
  index[orb_t::xy][orb_t::yz][orb_t::xy][orb_t::yz] = 40;
  index[orb_t::xz][orb_t::yz][orb_t::xz][orb_t::yz] = 40;
  index[orb_t::x2y2][orb_t::xz][orb_t::x2y2][orb_t::xz] = 40;
  index[orb_t::x2y2][orb_t::yz][orb_t::x2y2][orb_t::yz] = 40;
  expr[40] = OneCenterTwoElectronIntegralExpression(3.0 / 49.0, F2dd, 20.0 / 441.0, F4dd);

  index[orb_t::x2y2][orb_t::x2y2][orb_t::xy][orb_t::xy] = 41;
  expr[41] = OneCenterTwoElectronIntegralExpression(1, F0dd, 4.0 / 49.0, F2dd, -34.0 / 441.0, F4dd);

  index[orb_t::x2y2][orb_t::xz][orb_t::xy][orb_t::yz] = 42;
  expr[42] = OneCenterTwoElectronIntegralExpression(3.0 / 49.0, F2dd, -15.0 / 441.0, F4dd);

  index[orb_t::x2y2][orb_t::yz][orb_t::xy][orb_t::xz] = 43;
  expr[43] = OneCenterTwoElectronIntegralExpression(-3.0 / 49.0, F2dd, 15.0 / 441.0, F4dd);

  index[orb_t::xy][orb_t::xy][orb_t::xz][orb_t::xz] = 44;
  index[orb_t::xy][orb_t::xy][orb_t::yz][orb_t::yz] = 44;
  index[orb_t::xz][orb_t::xz][orb_t::yz][orb_t::yz] = 44;
  index[orb_t::x2y2][orb_t::x2y2][orb_t::xz][orb_t::xz] = 44;
  index[orb_t::x2y2][orb_t::x2y2][orb_t::yz][orb_t::yz] = 44;
  expr[44] = OneCenterTwoElectronIntegralExpression(1, F0dd, -2.0 / 49.0, F2dd, -4.0 / 441.0, F4dd);

  index[orb_t::s][orb_t::z2][orb_t::z2][orb_t::z2] = 45;
  expr[45] = OneCenterTwoElectronIntegralExpression(2.0 / sqrt(245), R2sddd);

  index[orb_t::s][orb_t::x2y2][orb_t::z2][orb_t::x2y2] = 46;
  index[orb_t::s][orb_t::xy][orb_t::z2][orb_t::xy] = 46;
  index[orb_t::s][orb_t::z2][orb_t::x2y2][orb_t::x2y2] = 46;
  index[orb_t::s][orb_t::z2][orb_t::xy][orb_t::xy] = 46;
  expr[46] = OneCenterTwoElectronIntegralExpression(-2.0 / sqrt(245), R2sddd);

  index[orb_t::s][orb_t::xz][orb_t::z2][orb_t::xz] = 47;
  index[orb_t::s][orb_t::yz][orb_t::z2][orb_t::yz] = 47;
  index[orb_t::s][orb_t::z2][orb_t::xz][orb_t::xz] = 47;
  index[orb_t::s][orb_t::z2][orb_t::yz][orb_t::yz] = 47;
  expr[47] = OneCenterTwoElectronIntegralExpression(1.0 / sqrt(245), R2sddd);

  index[orb_t::s][orb_t::xz][orb_t::x2y2][orb_t::xz] = 48;
  index[orb_t::s][orb_t::xz][orb_t::xy][orb_t::yz] = 48;
  index[orb_t::s][orb_t::yz][orb_t::xy][orb_t::xz] = 48;
  index[orb_t::s][orb_t::xy][orb_t::xz][orb_t::yz] = 48;
  index[orb_t::s][orb_t::x2y2][orb_t::xz][orb_t::xz] = 48;
  expr[48] = OneCenterTwoElectronIntegralExpression(sqrt(3) / sqrt(245), R2sddd);

  index[orb_t::s][orb_t::yz][orb_t::x2y2][orb_t::yz] = 49;
  index[orb_t::s][orb_t::x2y2][orb_t::yz][orb_t::yz] = 49;
  expr[49] = OneCenterTwoElectronIntegralExpression(-sqrt(3) / sqrt(245), R2sddd);

  index[orb_t::s][orb_t::z2][orb_t::x][orb_t::x] = 50;
  index[orb_t::s][orb_t::z2][orb_t::y][orb_t::y] = 50;
  expr[50] = OneCenterTwoElectronIntegralExpression(-1.0 / sqrt(125), R2sdpp);

  index[orb_t::s][orb_t::y][orb_t::y][orb_t::x2y2] = 51;
  expr[51] = OneCenterTwoElectronIntegralExpression(-sqrt(3) / sqrt(45), R1sppd);

  index[orb_t::s][orb_t::x2y2][orb_t::y][orb_t::y] = 52;
  expr[52] = OneCenterTwoElectronIntegralExpression(-sqrt(3) / sqrt(125), R2sdpp);

  index[orb_t::s][orb_t::z][orb_t::z][orb_t::z2] = 53;
  expr[53] = OneCenterTwoElectronIntegralExpression(2.0 / sqrt(45), R1sppd);

  index[orb_t::s][orb_t::z][orb_t::x][orb_t::xz] = 54;
  index[orb_t::s][orb_t::x][orb_t::z][orb_t::xz] = 54;
  index[orb_t::s][orb_t::z][orb_t::y][orb_t::yz] = 54;
  index[orb_t::s][orb_t::y][orb_t::z][orb_t::yz] = 54;
  index[orb_t::s][orb_t::x][orb_t::x][orb_t::x2y2] = 54;
  index[orb_t::s][orb_t::x][orb_t::y][orb_t::xy] = 54;
  index[orb_t::s][orb_t::y][orb_t::x][orb_t::xy] = 54;
  expr[54] = OneCenterTwoElectronIntegralExpression(sqrt(3) / sqrt(45), R1sppd);

  index[orb_t::s][orb_t::z2][orb_t::z][orb_t::z] = 55;
  expr[55] = OneCenterTwoElectronIntegralExpression(2.0 / sqrt(125), R2sdpp);

  index[orb_t::s][orb_t::xz][orb_t::x][orb_t::z] = 56;
  index[orb_t::s][orb_t::yz][orb_t::y][orb_t::z] = 56;
  index[orb_t::s][orb_t::x2y2][orb_t::x][orb_t::x] = 56;
  index[orb_t::s][orb_t::xy][orb_t::x][orb_t::y] = 56;
  expr[56] = OneCenterTwoElectronIntegralExpression(sqrt(3) / sqrt(125), R2sdpp);

  index[orb_t::s][orb_t::x][orb_t::x][orb_t::z2] = 57;
  index[orb_t::s][orb_t::y][orb_t::y][orb_t::z2] = 57;
  expr[57] = OneCenterTwoElectronIntegralExpression(-1.0 / sqrt(45), R1sppd);
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
