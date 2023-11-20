/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_OBSERVERS_H
#define SPARROW_OBSERVERS_H

#include <Utils/Typenames.h>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

class GradientObserver {
 public:
  void notifyGradient(const Utils::GradientCollection& gradient);

 private:
  virtual void notifyGradientImpl(const Utils::GradientCollection& gradient) = 0;
};

inline void GradientObserver::notifyGradient(const Utils::GradientCollection& gradient) {
  notifyGradientImpl(gradient);
}

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_OBSERVERS_H
