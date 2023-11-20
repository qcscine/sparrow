/**
 * @file ARealTimeSpectroscopyTest.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/RealTimeSpectroscopy/Utils/GuessPropagator.h>
#include <Utils/Typenames.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

class ARealTimeSpectroscopyTest : public Test {
 private:
  void SetUp() final {
  }
};

TEST_F(ARealTimeSpectroscopyTest, CanConstructGuessPropagator) {
  GuessPropagator guessProp;
}

TEST_F(ARealTimeSpectroscopyTest, CanRecordNewPosition) {
  Utils::PositionCollection pos = Utils::PositionCollection::Random(5, 3);
  GuessPropagator guessProp;
  guessProp.record(pos);
}

TEST_F(ARealTimeSpectroscopyTest, CanRecordNewSolution) {
  Eigen::MatrixXd lastSolutions = Eigen::MatrixXd::Random(100, 10);
  GuessPropagator guessProp;
  guessProp.record(lastSolutions);
}

TEST_F(ARealTimeSpectroscopyTest, CanRecordNewSolutionAndPositions) {
  Utils::PositionCollection pos = Utils::PositionCollection::Random(5, 3);
  Eigen::MatrixXd lastSolutions = Eigen::MatrixXd::Random(100, 10);
  GuessPropagator guessProp;
  guessProp.record(pos);
  guessProp.record(lastSolutions);
}

TEST_F(ARealTimeSpectroscopyTest, GetsEmptyVectorIfNothingRecorded) {
  Utils::PositionCollection pos = Utils::PositionCollection::Random(5, 3);
  GuessPropagator guessProp;
  Eigen::MatrixXd newGuess = guessProp.calculateGuessAtNewPosition(pos);
  ASSERT_EQ(newGuess.size(), 0);
}

TEST_F(ARealTimeSpectroscopyTest, GetsSameVectorForSequenceOfEqualStructures) {
  Utils::PositionCollection pos = Utils::PositionCollection::Random(5, 3);
  Eigen::MatrixXd lastSolutions = Eigen::MatrixXd::Random(100, 10);
  GuessPropagator guessProp;
  guessProp.record(pos);
  guessProp.record(lastSolutions);
  guessProp.record(pos);
  guessProp.record(lastSolutions);
  guessProp.record(pos);
  guessProp.record(lastSolutions);
  guessProp.record(pos);
  guessProp.record(lastSolutions);
  Eigen::MatrixXd newGuess = guessProp.calculateGuessAtNewPosition(pos);

  ASSERT_THAT(newGuess.rows(), lastSolutions.rows());
  ASSERT_THAT(newGuess.cols(), lastSolutions.cols());

  for (int row = 0; row < lastSolutions.rows(); ++row)
    for (int col = 0; col < lastSolutions.cols(); ++col)
      EXPECT_THAT(newGuess(row, col), DoubleNear(lastSolutions(row, col), 1.0e-12));
}

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine
