/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/VuvB.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ChargeSeparationParameter.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/KlopmanParameter.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;
using namespace multipole;

class AVuvBMatrix : public Test {
 public:
  AVuvBMatrix() : v(2) {
  }
  VuvB v;
  KlopmanParameter KV;
  ChargeSeparationParameter DV;
  const double pcore{3.923};
  const double zb{4.00};

  void SetUp() override {
    KV.set(MultipolePair::ss0, 2.27401429);
    KV.set(MultipolePair::sp1, 1.49936943);
    KV.set(MultipolePair::pp2, 1.99647866);
    KV.set(MultipolePair::sd2, 1.34719346);
    KV.set(MultipolePair::pd1, 1.88688065);
    KV.set(MultipolePair::dd2, 1.29691630);
    KV.set(MultipolePair::pp0, 2.27401429);
    KV.set(MultipolePair::dd0, 1.76187597);
    DV.computeFromExponents(4, 4, 3, 1.974330, 1.063106, 1.394806);
  }
};

TEST_F(AVuvBMatrix, ReturnsSameValuesAs2e2cMatrix) {
  Eigen::Vector3d Rab(0.0837, -0.6, 1.1);
  v.calculate<Utils::DerivativeOrder::Zero>(Rab, DV, KV, pcore, zb);

  KlopmanParameter KCore;
  ChargeSeparationParameter DCore;
  KCore.set(MultipolePair::ss0, pcore);
  Global2c2eMatrix g(2, 0, DV, DCore, KV, KCore);
  g.calculate<Utils::DerivativeOrder::Zero>(Rab);

  for (int i = 0; i < 40; i++) {
    SCOPED_TRACE("... for i = " + std::to_string(static_cast<long long>(i)));
    EXPECT_THAT(v.get(i), DoubleEq(-zb * g.get(i, 0)));
  }
}
} // namespace Sparrow
} // namespace Scine
