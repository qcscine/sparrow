/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_SPECTRUM_H
#define SPARROW_SPECTRUM_H

#include <Eigen/Core>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

class IncorrectSpectrumSizeException : public std::exception {
  const char* what() const noexcept final {
    return "Incorrect number of data points in spectrum assignment.";
  }
};

class Spectrum {
 public:
  /**
   * @brief Constructor initializing the number of data points in the spectrum.
   * @param size The number of (x,y) points in the spectrum.
   */
  explicit Spectrum(int size) {
    xData_.resize(size);
    yData_.resize(size);
  }

  /**
   * @brief Constructor initializing the x and y points.
   * @param xData
   * @param yData
   */
  Spectrum(const Eigen::VectorXd& xData, const Eigen::VectorXd& yData) {
    if (xData.size() != yData.size())
      throw IncorrectSpectrumSizeException();
    xData_ = xData;
    yData_ = yData;
  }
  Spectrum(const Eigen::VectorXd& xData, const Eigen::VectorXd& yData, std::pair<std::string, std::string> labels)
    : Spectrum(xData, yData) {
    labels_ = std::move(labels);
  }

  /**
   * @brief Getter for the number of (x,y) points in the spectrum.
   */
  int size() const {
    return xData_.size();
  }

  /**
   * @brief Gets the x data points in the spectrum.
   */
  const Eigen::VectorXd& getXData() const {
    return xData_;
  }

  /**
   * @brief Gets the y data points in the spectrum.
   */
  const Eigen::VectorXd& getYData() const {
    return yData_;
  }

  /**
   * @brief Gets the x data points in the spectrum.
   */
  double getXData(int point) const {
    return xData_(point);
  }

  /**
   * @brief Gets the y data points in the spectrum.
   */
  double getYData(int point) const {
    return yData_(point);
  }

  const std::pair<std::string, std::string>& getLabels() const {
    return labels_;
  }
  /**
   * @brief Changes the x data points in the spectrum.
   * @throws If xData does not have the right dimension.
   */
  void setXData(const Eigen::VectorXd& xData) {
    if (xData.size() != xData_.size())
      throw IncorrectSpectrumSizeException();
    xData_ = xData;
  }

  /**
   * @brief Changes the y data points in the spectrum.
   * @throws If yData does not have the right dimension.
   */
  void setYData(const Eigen::VectorXd& yData) {
    if (yData.size() != yData_.size())
      throw IncorrectSpectrumSizeException();
    yData_ = yData;
  }

  /**
   * @brief Changes the label of the x axis.
   */
  void setXLabel(std::string xLabel) {
    labels_.first = std::move(xLabel);
  }
  /**
   * @brief Changes the label of the y axis.
   */
  void setYLabel(std::string yLabel) {
    labels_.second = std::move(yLabel);
  }

 private:
  Eigen::VectorXd xData_;
  Eigen::VectorXd yData_;
  std::pair<std::string, std::string> labels_;
};

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_SPECTRUM_H
