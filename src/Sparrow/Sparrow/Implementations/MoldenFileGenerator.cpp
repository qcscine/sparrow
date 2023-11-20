#include <utility>

/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "GenericMethodWrapper.h"
#include "MoldenFileGenerator.h"
#include <Sparrow/Implementations/Nddo/Pm6/Wrapper/PM6MethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Utils/NDDOInitializer.h>
#include <Sparrow/Implementations/Sto6gParameters.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/IO/TurbomoleMinimalBasisfile.h>
#include <Utils/Scf/MethodInterfaces/LcaoMethod.h>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <iomanip>
#include <iostream>
namespace Scine {
namespace Sparrow {
using namespace nddo::GeneralTypes;

std::map<int, nddo::GeneralTypes::orb_t> MoldenFileGenerator::indexMap_ = {
    // s
    {0, orb_t::s},
    // px, py, pz
    {1, orb_t::x},
    {2, orb_t::y},
    {3, orb_t::z},
    // z2, xz, yz, x2-y2, xy
    {7, orb_t::x2y2},
    {5, orb_t::xz},
    {4, orb_t::z2},
    {6, orb_t::yz},
    {8, orb_t::xy}};

MoldenFileGenerator::MoldenFileGenerator(const GenericMethodWrapper& calculator) : calculator_(calculator) {
}

void MoldenFileGenerator::generateWavefunctionInformation(std::ostream& out) const {
  out << "[Molden Format]" << std::endl;
  out << "Written by Sparrow, the semiempirical library "
         "of the SCINE software package"
      << std::endl;
  out.imbue(std::locale("C"));
  generateAtomBlock(out);
  generateGTOBlock(out);
  generateMolecularOrbitalsBlock(out);
}

void MoldenFileGenerator::generateAtomBlock(std::ostream& out) const {
  auto structure = *calculator_.getStructure();
  const auto& elements = structure.getElements();
  const auto& positions = structure.getPositions();
  out << "[Atoms] (AU)" << std::fixed << std::setprecision(10) << std::endl;
  for (int atom = 0; atom < structure.size(); ++atom) {
    const auto e = elements[atom];
    const auto& p = positions.row(atom);
    out << std::left << std::setw(5) << Utils::ElementInfo::symbol(e) << std::setw(10) << atom + 1 << std::setw(10)
        << Utils::ElementInfo::Z(e) << std::right << " " << std::setw(18) << p.x() << " " << std::setw(18) << p.y()
        << " " << std::setw(18) << p.z() << std::endl;
  }
}

void MoldenFileGenerator::generateGTOBlock(std::ostream& out) const {
  std::unordered_map<int, Utils::AtomicGtos> expansionMap;
  if (calculator_.name() == "MNDO") {
    expansionMap = Sto6g::mndo();
  }
  else if (calculator_.name() == "AM1") {
    expansionMap = Sto6g::am1();
  }
  else if (calculator_.name() == "RM1") {
    expansionMap = Sto6g::rm1();
  }
  else if (calculator_.name() == "PM3") {
    expansionMap = Sto6g::pm3();
  }
  else if (calculator_.name() == "PM6") {
    expansionMap = Sto6g::pm6();
  }
  else if (calculator_.name().find("DFTB") != std::string::npos) {
    expansionMap = Sto6g::dftb();
  }
  else {
    auto turbomoleBasisfile = calculator_.getStoNGExpansionPath();
    expansionMap = Utils::readTurbomoleBasisfile(turbomoleBasisfile);
  }

  bool smallBasisMethod = calculator_.name() == "MNDO" || calculator_.name() == "AM1" || calculator_.name() == "PM3" ||
                          calculator_.name() == "RM1";
  out << std::right << "[GTO]" << std::endl;
  const auto& structurePtr = calculator_.getStructure();
  const int N = structurePtr->size();
  for (int atom = 0; atom < N; ++atom) {
    out << std::left << atom + 1 << " " << 0 << std::endl;
    const auto& expansion = expansionMap.at(Utils::ElementInfo::Z(structurePtr->getElement(atom)));
    if (expansion.s) {
      const auto& s = expansion.s.value();
      out << std::left << std::setw(5) << "s" << std::setw(5) << s.gtfs.size() << std::setw(5) << "1.00" << std::endl;
      for (const auto& primitive : s.gtfs) {
        out << std::right << std::setw(18) << primitive.exponent << " " << std::setw(18) << primitive.coefficient << std::endl;
      }
    }
    if (expansion.p) {
      const auto& p = expansion.p.value();
      out << std::left << std::setw(5) << "p" << std::setw(5) << p.gtfs.size() << std::setw(5) << "1.00" << std::endl;
      for (const auto& primitive : p.gtfs) {
        out << std::right << std::setw(18) << primitive.exponent << " " << std::setw(18) << primitive.coefficient << std::endl;
      }
    }
    if (expansion.d && !smallBasisMethod) {
      const auto& d = expansion.d.value();
      out << std::left << std::setw(5) << "d" << std::setw(5) << d.gtfs.size() << std::setw(5) << "1.00" << std::endl;
      for (const auto& primitive : d.gtfs) {
        out << std::right << std::setw(18) << primitive.exponent << " " << std::setw(18) << primitive.coefficient << std::endl;
      }
    }
    out << std::endl;
  }
}

void MoldenFileGenerator::generateMolecularOrbitalsBlock(std::ostream& out) const {
  const auto& method = calculator_.getLcaoMethod();
  const auto& moEnergies = method.getSingleParticleEnergies();
  const auto& occupation = method.getElectronicOccupation();
  out << "[5D]" << std::endl;
  out << "[7F]" << std::endl;
  out << "[9G]" << std::endl;
  out << "[MO]" << std::endl;
  if (method.unrestrictedCalculationRunning()) {
    writeMOBlock(out, method.getMolecularOrbitals().alphaMatrix(), occupation.getFilledAlphaOrbitals(),
                 moEnergies.getAlphaEnergies(), "Alpha");
    writeMOBlock(out, method.getMolecularOrbitals().betaMatrix(), occupation.getFilledBetaOrbitals(),
                 moEnergies.getBetaEnergies(), "Beta");
  }
  else {
    writeMOBlock(out, method.getMolecularOrbitals().restrictedMatrix(), occupation.getFilledRestrictedOrbitals(),
                 moEnergies.getRestrictedEnergies(), "Alpha");
  }
}

void MoldenFileGenerator::writeMOBlock(std::ostream& out, Eigen::MatrixXd moMatrix, const std::vector<int>& filledOrb,
                                       const std::vector<double>& moEnergies, const std::string& spin) const {
  const auto& method = calculator_.getLcaoMethod();
  std::string occupancy = method.unrestrictedCalculationRunning() ? "1.0" : "2.0";
  for (int orb = 0; orb < moMatrix.cols(); ++orb) {
    const auto& mo = moMatrix.col(orb);
    bool isOccupied = std::find(filledOrb.begin(), filledOrb.end(), orb) != filledOrb.end();
    out << std::left << std::setw(10) << "Sym=" << std::setw(10) << "A1" << std::endl;
    out << std::left << std::setw(10) << "Ene=" << std::setw(10) << moEnergies[orb] << std::endl;
    out << std::left << std::setw(10) << "Spin=" << std::setw(10) << spin << std::endl;
    out << std::left << std::setw(10) << "Occup=" << std::setw(10) << (isOccupied ? occupancy : "0.0") << std::endl;
    for (int atom = 0; atom < method.getAtomsOrbitalsIndexesHolder().getNAtoms(); ++atom) {
      int firstAtom = method.getAtomsOrbitalsIndexesHolder().getFirstOrbitalIndex(atom);
      int nAOs = method.getAtomsOrbitalsIndexesHolder().getNOrbitals(atom);
      for (int ao = 0; ao < nAOs; ++ao) {
        out << std::right << std::setw(4) << firstAtom + ao + 1 << std::setw(18)
            << mo(firstAtom + indexMap_[static_cast<orb_t>(ao)]) << std::endl;
      }
    }
  }
}

} // namespace Sparrow
} // namespace Scine
