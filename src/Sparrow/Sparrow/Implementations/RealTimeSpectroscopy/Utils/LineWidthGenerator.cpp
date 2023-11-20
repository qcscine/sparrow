/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#define _USE_MATH_DEFINES
#include "LineWidthGenerator.h"
#include <Sparrow/Implementations/RealTimeSpectroscopy/Utils/Spectrum.h>
#include <algorithm>
#include <cmath>
#include <vector>
namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

LineWidthGenerator::LineWidthGenerator(const Spectrum& inputSpectrum) : inputSpectrum_(inputSpectrum) {
  assert(inputSpectrum_.size() == inputSpectrum.size());
}

auto LineWidthGenerator::lorentzFunction(double wavenumber, double intensity, double fwhm) const {
  assert(fwhm != 0 && "Lorentzian cannot be calculated if full width at half maximum is equal to 0!");
  return [=](double x) -> double {
    double offset = x - wavenumber;
    double lorentzDenominator = offset * offset + fwhm * fwhm * 0.25;
    double lorentzNumerator = 0.5 * fwhm / M_PI;
    return intensity * lorentzNumerator / lorentzDenominator;
  };
}

auto LineWidthGenerator::barplotFunction(const Eigen::VectorXd& waveNumbersXData, int normalModeIndex, double resolution) const {
  return [&](int point) -> double {
    auto peakPosition = inputSpectrum_.getXData(normalModeIndex) - waveNumbersXData[0];
    auto pointNumberOnVector = peakPosition / resolution;
    double intensity = static_cast<int>(pointNumberOnVector) == point ? inputSpectrum_.getYData(normalModeIndex) : 0.0;
    return intensity;
  };
}

Spectrum LineWidthGenerator::generateLorentzianProfile(double resolution, double fwhm) const {
  if (inputSpectrum_.getYData().size() == 0) {
    throw(std::runtime_error("Spectrum with no data!"));
  }
  Eigen::VectorXd wavenumberXVector = generateWavenumberXVector(resolution, fwhm);
  Spectrum spectrum(wavenumberXVector.size());
  Eigen::VectorXd intensityYVector = Eigen::VectorXd::Zero(wavenumberXVector.size());
  spectrum.setXData(wavenumberXVector);

  for (int normalModeIndex = 0; normalModeIndex < inputSpectrum_.getXData().size(); ++normalModeIndex) {
    // exclude roto/translational modes in case of IR spectrum
    if (inputSpectrum_.getXData(normalModeIndex) != 0) {
      std::function<double(double)> convolution =
          lorentzFunction(inputSpectrum_.getXData(normalModeIndex), inputSpectrum_.getYData(normalModeIndex), fwhm);
      if (fwhm == 0.)
        convolution = barplotFunction(wavenumberXVector, normalModeIndex, resolution);

      for (int point = 0; point < wavenumberXVector.size(); ++point) {
        intensityYVector(point) += convolution(wavenumberXVector(point));
      }
    }
  }

  spectrum.setYData(intensityYVector);
  spectrum.setXLabel(inputSpectrum_.getLabels().first);
  spectrum.setYLabel(inputSpectrum_.getLabels().second);
  return spectrum;
}

Eigen::VectorXd LineWidthGenerator::generateWavenumberXVector(double resolution, double fwhm) const {
  // Set first point to be sufficiently far to show the first peak
  double buffer = 3 * fwhm;
  std::pair<double, double> bounds = {inputSpectrum_.getXData().minCoeff() - buffer,
                                      inputSpectrum_.getXData().maxCoeff() + buffer};
  // Round up the number of points to avoid cropping in case of a line plot (FWHM = 0)
  auto numberOfPointOnXAxis = static_cast<unsigned>(std::ceil((bounds.second - bounds.first) / resolution));

  Eigen::VectorXd wavenumberXAxis = Eigen::VectorXd::Zero(numberOfPointOnXAxis);
  for (unsigned i = 0; i < wavenumberXAxis.size(); ++i) {
    wavenumberXAxis(i) = bounds.first + i * resolution;
  }
  return wavenumberXAxis;
}

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine
