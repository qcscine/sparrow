/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef RTSPECTROSCOPY_LINEWIDTHGENERATOR_H
#define RTSPECTROSCOPY_LINEWIDTHGENERATOR_H

#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

class Spectrum;

class LineWidthGenerator {
 public:
  explicit LineWidthGenerator(const Spectrum& inputData);

  Spectrum generateLorentzianProfile(double resolution, double fwhm) const;

 private:
  auto barplotFunction(const Eigen::VectorXd& waveNumbers, int normalModeIndex, double resolution) const;
  auto lorentzFunction(double wavenumber, double intensity, double fwhm) const;
  Eigen::VectorXd generateWavenumberXVector(double resolution, double fwhm) const;

  const Spectrum& inputSpectrum_;
};

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine

#endif // RTSPECTROSCOPY_LINEWIDTHGENERATOR_H
