#include <utility>

/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "GenericMethodWrapper.h"
#include "MoldenFileGenerator.h"
#include <Sparrow/Implementations/Nddo/Pm6/Wrapper/PM6MethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Utils/NDDOInitializer.h>
#include <Utils/IO/NativeFilenames.h>
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
  generateAtomBlock(out);
  generateGTOBlock(out);
  generateMolecularOrbitalsBlock(out);
}

std::vector<Utils::AtomicGtos> MoldenFileGenerator::getStoNGExpansion() const {
  auto filename = calculator_.getStoNGExpansionPath();
  auto elements = calculator_.getStructure()->getElements();
  std::vector<Utils::AtomicGtos> atomicGtosVector;

  for (Utils::ElementType el : elements) {
    std::ifstream basisFile(filename);
    if (!basisFile.is_open()) {
      throw std::runtime_error("Not able to open basis set file '" + filename + "'.");
    }
    // TODO: parse basis set file and give it back as is.
    std::string line;
    while (std::getline(basisFile, line)) {
      std::string symbol = Utils::ElementInfo::symbol(el);
      std::transform(line.begin(), line.end(), line.begin(), ::tolower);
      std::transform(symbol.begin(), symbol.end(), symbol.begin(), ::tolower);
      if (line.substr(0, symbol.size() + 1) == symbol + " ") {
        Utils::AtomicGtos atomicGtos;
        // found element
        std::string gtoLine;
        std::stringstream gtoStream;
        // skip asterisk
        std::getline(basisFile, gtoLine);

        while (std::getline(basisFile, gtoLine)) {
          if (gtoLine.find('*') != std::string::npos)
            break;
          // get first line
          int numberOfPrimitives = 0;
          std::string label;
          // read GTO contraction infos
          gtoStream = std::stringstream(gtoLine);
          gtoStream >> numberOfPrimitives >> label;

          int angularMomentum;
          if (label == "s")
            angularMomentum = 0;
          else if (label == "p")
            angularMomentum = 1;
          else if (label == "d")
            angularMomentum = 2;
          else
            throw std::runtime_error("Angular momentum not available in semiempirical methods.");
          Utils::GtoExpansion gto(numberOfPrimitives);
          gto.setAngularMomentum(angularMomentum);
          for (int primitive = 0; primitive < numberOfPrimitives; ++primitive) {
            double coefficient = 0.0, exponent = 0.0;
            std::getline(basisFile, gtoLine);
            gtoStream = std::stringstream(gtoLine);
            gtoStream >> exponent >> coefficient;
            gto.setGTF(primitive, exponent, coefficient);
          }
          if (angularMomentum == 0) {
            atomicGtos.setS(std::move(gto));
          }
          else if (angularMomentum == 1) {
            atomicGtos.setP(std::move(gto));
          }
          else {
            atomicGtos.setD(std::move(gto));
          }
        }
        atomicGtosVector.emplace_back(std::move(atomicGtos));
      }
    }
  }

  return atomicGtosVector;
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
  const auto& gtoExpansion = getStoNGExpansion();
  auto smallBasisMethod = calculator_.name() == "MNDO" || calculator_.name() == "AM1" || calculator_.name() == "PM3" ||
                          calculator_.name() == "RM1";
  out << std::right << "[GTO]" << std::endl;
  for (int atom = 0; atom < calculator_.getStructure()->size(); ++atom) {
    out << std::left << atom + 1 << " " << 0 << std::endl;
    if (gtoExpansion[atom].hasS()) {
      const auto& s = gtoExpansion[atom].s();
      out << std::left << std::setw(5) << "s" << std::setw(5) << s.size() << std::setw(5) << "1.00" << std::endl;
      for (int primitive = 0; primitive < s.size(); ++primitive) {
        out << std::right << std::setw(18) << s.getExponent(primitive) << " " << std::setw(18)
            << s.getCoefficient(primitive) << std::endl;
      }
    }
    if (gtoExpansion[atom].hasP()) {
      const auto& p = gtoExpansion[atom].p();
      out << std::left << std::setw(5) << "p" << std::setw(5) << p.size() << std::setw(5) << "1.00" << std::endl;
      for (int primitive = 0; primitive < p.size(); ++primitive) {
        out << std::right << std::setw(18) << p.getExponent(primitive) << " " << std::setw(18)
            << p.getCoefficient(primitive) << std::endl;
      }
    }
    if (gtoExpansion[atom].hasD() && !smallBasisMethod) {
      const auto& d = gtoExpansion[atom].d();
      out << std::left << std::setw(5) << "d" << std::setw(5) << d.size() << std::setw(5) << "1.00" << std::endl;
      for (int primitive = 0; primitive < d.size(); ++primitive) {
        out << std::right << std::setw(18) << d.getExponent(primitive) << " " << std::setw(18)
            << d.getCoefficient(primitive) << std::endl;
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
