/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/MolecularTrajectory.h>

int main() {
  auto traj =
      Scine::Utils::MolecularTrajectoryIO::read(Scine::Utils::MolecularTrajectoryIO::format::xyz,
                                                "/scratch/severinp/workdir/Report/BergmanTraj/trajBergmanAll.xyz");
  auto el = traj.getElementTypes();
  int i = 0;
  for (auto struc : traj) {
    std::string filename = "/scratch/severinp/workdir/Report/BergmanTraj/trajBergmanAll" + std::to_string(i) + ".xyz";
    Scine::Utils::ChemicalFileHandler::write(filename, {el, struc});
    ++i;
  }
  return 0;
}
