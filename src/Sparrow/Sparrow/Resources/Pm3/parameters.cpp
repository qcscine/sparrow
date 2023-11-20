/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 *
 * @note This file was generated by the Embed binary from runtime values.
 *   Prefer improving the generator over editing this file whenever possible.
 *
 * This file contains functions generating runtime values. It was
 * generated from runtime values of its return type. It is not intended to be
 * human-readable. A small guide: Return values are directly brace-initialized
 * in deep nesting to keep file size to a minimum. Types are annotated only when
 * necessary. Floating point values are represented in hexadecimal (see
 * std::hexfloat) to ensure serialization does not cause loss of accuracy.
 *
 * The functions defined here might be declared and called elsewhere entirely.
 *
 */
#include "Sparrow/Implementations/Nddo/Parameters.h"

namespace Scine {
namespace Sparrow {
namespace nddo {

Parameters pm3() {
  // clang-format off
  using Pack = Parameters::Atomic::Pack;
  using Spd = Parameters::Atomic::Pack::Spd;
  using Gauss = Parameters::Atomic::GaussianRepulsion;

  return {
    {
      {32, {Pack {Spd {-0x1.1bbcd0fe8ab5p+5, -0x1.f961b93da3cf8p+4, 0x0p+0}, Spd {-0x1.54ccd6ddc7c6dp+2, -0x1.2005227eb0094p+1, 0x0p+0}, Spd {0x1.1e6191eb4f6edp+1, 0x1.97a99df39b6dbp+0, 0x0p+0}, Spd {}, 0x1.58202b841248dp+2, 0x1.eaffd4cc758fp+2, 0x1.46b476ca61b88p+3, 0x1.bb272dd6d2e01p+2, 0x1.5646f80c15b8p+0, 0x1.43e311f31c97cp+1, 0x0p+0, 0x0p+0, 0x1.f8eb13dfb0d52p+0}, {{Gauss {0x1.ed24f582ce7cp-1, 0x1.80c4d3da0747fp+2, 0x1.14e92923e5b85p+1}, Gauss {-0x1.eb350c51471a7p-1, 0x1.6ff2918273e32p+2, 0x1.15adfeab6c29fp+1}}}}},
      {31, {Pack {Spd {-0x1.ddb082491afcp+4, -0x1.5e0185058dde8p+4, 0x0p+0}, Spd {-0x1.3c85015c20924p+2, -0x1.a0d2806af46aap-2, 0x0p+0}, Spd {0x1.d8d79d0a67621p+0, 0x1.adc74751ce28fp-1, 0x0p+0}, Spd {}, 0x1.0eac79702e664p+3, 0x1.458f08461f9fp+2, 0x1.1d9eabbcb1cc9p+3, 0x1.3eea35935fc3bp+2, 0x1.068fb00bcbe62p+1, 0x1.9bc77d579999fp+0, 0x0p+0, 0x0p+0, 0x1.9ae8d10f51acap+0}, {{Gauss {-0x1.1ecfc829cfdd2p-1, 0x1.67e3b46fdeb53p+2, 0x1.8822bbecaab8ap+0}, Gauss {-0x1.1746cb966be7bp-2, 0x1.fde96c3fc43b3p+0, 0x1.1788db0574b4p+1}}}}},
      {30, {Pack {Spd {-0x1.2883e20ccff22p+4, -0x1.61845fe111277p+3, 0x0p+0}, Spd {-0x1.6e603d577963ap-1, -0x1.9684f09528f19p+2, 0x0p+0}, Spd {0x1.d1eacc92146a2p+0, 0x1.81c5a3e39f773p+0, 0x0p+0}, Spd {}, 0x1.35ab96f21f6cbp+3, 0x1.3ebb2bba98edap+2, 0x1.ef1df761cbccfp+2, 0x1.2adba4d6e47dcp+2, 0x1.3369552e2fbe3p-1, 0x1.67ec98c5085d4p+0, 0x0p+0, 0x0p+0, 0x1.59a1db877ab32p+0}, {{Gauss {-0x1.c79d4d834091cp-4, 0x1.8018372e6a769p+2, 0x1.841aac53b0814p+0}, Gauss {-0x1.0f1800a7c5ac4p-3, 0x1.feef4e0114d2fp+0, 0x1.4283a109d0636p+1}}}}},
      {3, {Pack {Spd {-0x1.5333333333333p+2, -0x1.b333333333333p+1, 0x0p+0}, Spd {-0x1.199999999999ap-1, -0x1.8p+0, 0x0p+0}, Spd {0x1.4cccccccccccdp-1, 0x1.8p-1, 0x0p+0}, Spd {}, 0x1.2p+2, 0x1.5p+2, 0x1.8p+1, 0x1.2p+2, 0x1.3333333333333p-3, 0x1.83019f9257052p+1, 0x0p+0, 0x0p+0, 0x1.4147ae147ae14p+0}, {{Gauss {-0x1.ccccccccccccdp-2, 0x1.4p+2, 0x1p+0}, Gauss {0x1.999999999999ap-1, 0x1.ap+2, 0x1p+0}}}}},
      {50, {Pack {Spd {-0x1.1466cb10342abp+5, -0x1.9e4f8a4c1ebc8p+4, 0x0p+0}, Spd {-0x1.649528f190d17p+1, -0x1.00c493426783ap+1, 0x0p+0}, Spd {0x1.2fc9363f572dep+1, 0x1.a3633ce63a5c2p+0, 0x0p+0}, Spd {}, 0x1.4614c0160525p+3, 0x1.6b1fb3fa6defcp+2, 0x1.cf0f98fa37692p+2, 0x1.4ba964e8b7e4ep+2, 0x1.087cfa26a22b4p+0, 0x1.55cf685d8c782p+0, 0x0p+0, 0x0p+0, 0x1.b31c432ca57a8p+0}, {{Gauss {-0x1.33ec460ed80a1p-3, 0x1.805d4a5df223ap+2, 0x1.b4636b0963561p+0}, Gauss {-0x1.6bdd334c5da6ap-5, 0x1.20f1dc50ce4ebp+1, 0x1.3c24aada33bdap+1}}}}},
      {20, {Pack {Spd {-0x1.7ebf005b6f926p+3, -0x1.14529faa54d2cp+3, 0x0p+0}, Spd {-0x1.f94166d29cf06p-1, -0x1.fbf567f6895ecp+1, 0x0p+0}, Spd {0x1.3570153bd1676p+0, 0x1.e1c27e9531551p-1, 0x0p+0}, Spd {}, 0x1.7e3519bd4067dp+2, 0x1.dc59c644f34f1p+1, 0x1.3d7d63022dd7bp+2, 0x1.db17580b59d05p+1, 0x1.95ebd9018e758p-1, 0x1.239db6b4b279cp+1, 0x0p+0, 0x0p+0, 0x1.ffff73124653ap-2}, {{}}}},
      {2, {Pack {Spd {-0x1.1bf1e75360d02p+5, 0x1.3ffeb76053fdbp+3, 0x0p+0}, Spd {-0x1.18df3d7d3910cp+4, -0x1.a0d779cef8f6fp+4, 0x0p+0}, Spd {0x1.c5653e483a2b7p+0, 0x1.b9b7838f6c193p+2, 0x0p+0}, Spd {}, 0x1.01e3704790b85p+3, 0x1.36263b4db3e8p+4, 0x1.607cb4bc34cf1p+3, 0x1.6b414f55a8a9p+3, 0x1.107219220ff54p-1, 0x1.b031a843ab7b7p+0, 0x0p+0, 0x0p+0, 0x1.72025834747bdp+2}, {{Gauss {0x1.96ba56a883417p-3, 0x1.003ede6220869p+0, 0x1.019cfea8e7eb3p-1}}}}},
      {49, {Pack {Spd {-0x1.a2d1bc5586445p+4, -0x1.4017d8cf398e9p+4, 0x0p+0}, Spd {-0x1.7f2513b5bf6a1p+1, -0x1.d433508f648c7p+0, 0x0p+0}, Spd {0x1.021016ce789e7p+1, 0x1.72027525460aap+0, 0x0p+0}, Spd {}, 0x1.a3837b4a2339cp+2, 0x1.93273929ed395p+2, 0x1.075b1e9f27781p+3, 0x1.3efd50225742ep+2, 0x1.50d3b6cbd987cp+1, 0x1.09aef779c99a9p+1, 0x0p+0, 0x0p+0, 0x1.6b1b4784230fdp+0}, {{Gauss {-0x1.5f5f91600f345p-2, 0x1.fe7903211cb04p+0, 0x1.a021d10b1feebp+0}, Gauss {-0x1.c0a4a05dd8f93p-4, 0x1.6bb9d3cbc48f1p+2, 0x1.6efa26a22b389p+1}}}}},
      {36, {Pack {Spd {0x1.398ee61720adfp+3, -0x1.2770321881004p+6, 0x0p+0}, Spd {-0x1.359fdf2df860dp+0, -0x1.0dd9b837ea522p+3, 0x0p+0}, Spd {0x1.c7ca55234803ap+1, 0x1.fbb366ca39749p+0, 0x0p+0}, Spd {}, 0x1.ff6724d247601p-2, 0x1.3ffacad729472p+4, 0x1.30d605c6129e7p+3, 0x1.665af94ada18ep+3, 0x1.1746547827591p+1, 0x1.b3e3f5d63078cp+4, 0x0p+0, 0x0p+0, 0x1.b3d9e40c8472cp+0}, {{Gauss {0x1.8f95b35b0bbf5p-1, 0x1.1ac9d93235dd3p+3, 0x1.00057dee42686p-1}}}}},
      {83, {Pack {Spd {-0x1.0bf7ae5796bfdp+5, -0x1.1c2b0fadf2ecfp+5, 0x0p+0}, Spd {-0x1.66ddb9841aac5p+2, -0x1.7335b0bbf50e3p+2, 0x0p+0}, Spd {0x1.3aa7221858bc6p+2, 0x1.ef57e670e2c13p+0, 0x0p+0}, Spd {}, 0x1.3f53a3ec02f3p+2, 0x1.1645b078d92fbp+3, 0x1.869c99285a922p+2, 0x1.0abbfb58d1527p+3, 0x1.32c01e68a0d35p-1, 0x1.5d0a433a675ap+1, 0x0p+0, 0x0p+0, 0x1.db809917939a8p+0}, {{Gauss {0x1.4a74ea8da7f3dp+1, 0x1.460474d9c6b05p+2, 0x1.ffc829cfdd228p-2}, Gauss {0x1.ee2435696e58ap-5, 0x1.801932d6ece14p+2, 0x1.36c7b890d5a5cp+1}}}}},
      {13, {Pack {Spd {-0x1.8d86c6583e857p+4, -0x1.6439fec99f1aep+4, 0x0p+0}, Spd {-0x1.3048387df5cf2p-1, -0x1.e9c0ebedfa44p-1, 0x0p+0}, Spd {0x1.b3f077ccc0379p+0, 0x1.12d959a30984ep+0, 0x0p+0}, Spd {}, 0x1.71b60f1b25f63p+2, 0x1.964230fcf80dcp+2, 0x1.751d8a5482385p+3, 0x1.87bfb9bed30fp+2, 0x1.0066516db0dd8p+2, 0x1.2d79034c0fcep+1, 0x0p+0, 0x0p+0, 0x1.858e53eb399f6p+0}, {{Gauss {-0x1.e471b4784231p-2, 0x1.ea7381d7dbf48p+0, 0x1.73a4723aafff3p+0}, Gauss {-0x1.3b7f1737542a2p-3, 0x1.8053543aeab7ap+2, 0x1.428f42fe82518p+1}}}}},
      {35, {Pack {Spd {-0x1.d27a2ca9ac365p+6, -0x1.28e894812be49p+6, 0x0p+0}, Spd {-0x1.f2bdd11be6e65p+4, -0x1.b418c9fb6134dp+2, 0x0p+0}, Spd {0x1.564d1e96c3fc4p+2, 0x1.1054de7ea5f85p+1, 0x0p+0}, Spd {}, 0x1.fe3089a027525p+3, 0x1.090c64fdb09a6p+3, 0x1.00fca42aed139p+4, 0x1.f4474107314cbp+2, 0x1.2861847f56217p-1, 0x1.b4ed5558770f8p-1, 0x0p+0, 0x0p+0, 0x1.418409e55c0fdp+1}, {{Gauss {0x1.ebc126a65cf68p-1, 0x1.7e7f1b6912125p+2, 0x1.292bf5515054bp+1}, Gauss {-0x1.e8eabffcdab19p-1, 0x1.7c7603925bb7bp+2, 0x1.2a008e9b38d61p+1}}}}},
      {82, {Pack {Spd {-0x1.e52a023209678p+4, -0x1.86d0374ff865dp+4, 0x0p+0}, Spd {-0x1.8810c6f7a0b5fp+2, -0x1.653ae685db76bp+0, 0x0p+0}, Spd {0x1.9215c209246bfp+1, 0x1.e475818c5c9a3p+0, 0x0p+0}, Spd {}, 0x1.c0c47a17f4129p+2, 0x1.4bc30d306a2b1p+2, 0x1.b2cd530489d28p+2, 0x1.42ebf22c01e69p+2, 0x1.90f92af9a8cdfp+0, 0x1.f0ba86792b651p+0, 0x0p+0, 0x0p+0, 0x1.9ebb44e50c5ebp+0}, {{Gauss {-0x1.f61240746455fp-4, 0x1.80322af5771p+2, 0x1.e6cf0f9d2bf55p+0}, Gauss {-0x1.d00f776c4827bp-5, 0x1.2f98dcdb37c9ap+2, 0x1.6e520d130df9cp+1}}}}},
      {12, {Pack {Spd {-0x1.d3f540895d0b7p+3, -0x1.c58cfbfc6540dp+3, 0x0p+0}, Spd {-0x1.092d2bb23571dp+1, -0x1.23a01eeed8905p-1, 0x0p+0}, Spd {0x1.65a89b951c5c5p-1, 0x1.7bc393682730cp+0, 0x0p+0}, Spd {}, 0x1.ac6f694467382p+2, 0x1.ba44bf4cb1898p+2, 0x1.b2d0d0678c005p+2, 0x1.c5d00b45ae6p+2, 0x1.162b6ae7d566dp-1, 0x1.0426a5055d8dp+1, 0x0p+0, 0x0p+0, 0x1.5442fa5093965p+0}, {{Gauss {0x1.0efb7e90ff972p+1, 0x1.809b456b441bcp+2, 0x1.0acdd0d8cb07dp+1}, Gauss {-0x1.461d3aa369fcfp+1, 0x1.194dbdf8f473p+2, 0x1.08267839cd812p+1}}}}},
      {1, {Pack {Spd {-0x1.a258a54823854p+3, 0x0p+0, 0x0p+0}, Spd {-0x1.6818c5c9a34cap+2, 0x0p+0, 0x0p+0}, Spd {0x1.ef84662bae03bp-1, 0x0p+0, 0x0p+0}, Spd {}, 0x1.d96a26e547171p+3, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.d6de1f04a8057p-1, 0x0p+0, 0x0p+0, 0x1.ad9e0e736049fp+1}, {{Gauss {0x1.20f5c28f5c28fp+0, 0x1.46297bfa4c61ep+2, 0x1.89974e65bea0cp+0}, Gauss {-0x1.0f71b8aa00193p+0, 0x1.803e10060781p+2, 0x1.91f7e80389f84p+0}}}}},
      {48, {Pack {Spd {-0x1.fa83c297bfa4cp+3, 0x1.17fe52157689dp+3, 0x0p+0}, Spd {-0x1.129f4906034f4p+3, -0x1.33baba7b9170dp-1, 0x0p+0}, Spd {0x1.ade9f2778140ep+0, 0x1.0880303c07ee1p+1, 0x0p+0}, Spd {}, 0x1.269f6a93f290bp+3, 0x1.3cadbc664d3bfp+2, 0x1.0768c47a17f41p+3, 0x1.2adba4d6e47dcp+2, 0x1.a7fef39085f4ap+0, 0x1.7a4e950e6f847p+0, 0x0p+0, 0x0p+0, 0x1.867f6f4be835ep+0}, {{}}}},
      {34, {Pack {Spd {-0x1.bb066ba493c8ap+5, -0x1.8e95a8deb0faep+5, 0x0p+0}, Spd {-0x1.8a19c17225b75p+2, -0x1.5f8df37329c34p+2, 0x0p+0}, Spd {0x1.69fd933e35c5bp+1, 0x1.bb877ab324852p+0, 0x0p+0}, Spd {}, 0x1.dbaf922962cfep+2, 0x1.322fba01eeed9p+3, 0x1.41ef4be835deep+3, 0x1.ee5ac03ff6901p+2, 0x1.010f49491f2dcp+2, 0x1.d49e97854249ap+0, 0x0p+0, 0x0p+0, 0x1.85a0620ab7132p+1}, {{Gauss {0x1.882cf52b90a78p-5, 0x1.80793dd97f62bp+2, 0x1.0a75b3e1437c5p+1}, Gauss {0x1.d5e4a38327675p-4, 0x1.808e15011904bp+2, 0x1.84344c37e6f72p+0}}}}},
      {81, {Pack {Spd {-0x1.e0d9c8c9320dap+4, -0x1.aebaeddce7cdp+4, 0x0p+0}, Spd {-0x1.15a176ddaceeep+0, -0x1.fc985ad538ac2p+2, 0x0p+0}, Spd {0x1.b78c0485a0be5p+2, 0x1.f82d8c2a454dep+0, 0x0p+0}, Spd {}, 0x1.4ebbb1f255f35p+3, 0x1.3f89ca18bd662p+2, 0x1.672a0cae642cp+3, 0x1.1ecea8da7f3cfp+3, 0x1.43e45803cd142p+1, 0x1.4cf9a242e07a5p+0, 0x0p+0, 0x0p+0, 0x1.5748909289daep+0}, {{Gauss {-0x1.5c84a515ce9e6p+0, 0x1.c7532e7b3d8ep+1, 0x1.17c1df3300de5p+0}, Gauss {-0x1.73eccc46950fcp-5, 0x1.274b9cb6848bfp+1, 0x1.7b8611fd5885dp+1}}}}},
      {11, {Pack {Spd {-0x1.40e65965cce37p+2, -0x1.766df09e632cfp+1, 0x0p+0}, Spd {-0x1.f96cd799af691p+1, -0x1.0f5f1fb5a7ed2p+2, 0x0p+0}, Spd {0x1.54b8efa0366bdp+1, 0x1.c479e59f2ba9dp-1, 0x0p+0}, Spd {}, 0x1.730aa50a10244p+2, 0x1.398bdb6af5416p+0, 0x1.40000d009982fp+3, 0x1.ffd215e1e9948p-1, 0x1.9988e6c3eed7cp-2, 0x1.2c64492ac149bp+1, 0x0p+0, 0x0p+0, 0x1.d263ce05b1f0bp-1}, {{}}}},
      {5, {Pack {Spd {-0x1.93d24b698ade1p+5, -0x1.2b4bbe0157eedp+5, 0x0p+0}, Spd {-0x1.51975b9c0808ep+3, -0x1.fff2bd215e1eap+1, 0x0p+0}, Spd {0x1.8800a2bd2eca1p+0, 0x1.24b9c65fcb41bp+0, 0x0p+0}, Spd {}, 0x1.2473d54f524dbp+4, 0x1.8a1b82a7e58cbp+3, 0x1.eaa87cc11bbeap+3, 0x1.65b68f3df604dp+3, 0x1.33177a7008a69p-1, 0x1.7d1d44338b9dbp-1, 0x0p+0, 0x0p+0, 0x1.1aeeebdb85cd3p+1}, {{Gauss {-0x1.6848edaf9b63ap-2, 0x1.801c3fd1a7272p+1, 0x1.a5f2bdf81db37p-1}}}}},
      {52, {Pack {Spd {-0x1.67811904b3c3ep+5, -0x1.7283465625a68p+5, 0x0p+0}, Spd {-0x1.5523810e8859p+1, -0x1.f29d7342edbb6p+1, 0x0p+0}, Spd {0x1.0a976bc1effap+2, 0x1.a5c62a1b5c7cep+0, 0x0p+0}, Spd {}, 0x1.48298eda22f6ap+3, 0x1.f1c4113c68662p+2, 0x1.0569a2c669058p+3, 0x1.f053e707e175dp+2, 0x1.e2e008e9b38d6p+1, 0x1.53a470e317111p+0, 0x0p+0, 0x0p+0, 0x1.3e151a437824dp+1}, {{Gauss {0x1.118a009f62307p-5, 0x1.7d355043e5322p+2, 0x1.238793dd97f63p+1}, Gauss {-0x1.ebff79c842fa5p+0, 0x1.3e4938583622p+2, 0x1.0c6994185058ep-1}}}}},
      {10, {Pack {Spd {0x1.34f9a3c56c864p+3, -0x1.1c962abc6c2c8p+6, 0x0p+0}, Spd {-0x1.35caca940d68ap-3, -0x1.7d8d637af99b6p+4, 0x0p+0}, Spd {0x1.7ffdf19d66adbp+2, 0x1.0b3775b813016p+2, 0x0p+0}, Spd {}, 0x1.fee9721f014dp-2, 0x1.2f1eb5cdacc6ap+4, 0x1.451833d36c234p+3, 0x1.121352f8cb586p+3, 0x1.3332f6cd51571p-2, 0x1.b44f3376fade3p+4, 0x0p+0, 0x0p+0, 0x1.44d4306e5cd4fp+1}, {{Gauss {0x1.e1b640994d438p-3, 0x1.19f91fdc3e59dp+3, 0x1.17c0fb0772bb1p+0}}}}},
      {6, {Pack {Spd {-0x1.7a299d883ba34p+5, -0x1.2222a5e785b5bp+5, 0x0p+0}, Spd {-0x1.7d1ed7c6fbd27p+3, -0x1.39b02b40f66a5p+3, 0x0p+0}, Spd {0x1.90a9691a75cd1p+0, 0x1.d7a3ec02f2f98p+0, 0x0p+0}, Spd {}, 0x1.666c332f01755p+3, 0x1.597b395c42203p+3, 0x1.487b19e731d2ep+3, 0x1.215cb35f3d7d4p+3, 0x1.253ed527e5215p+1, 0x1.36f7b35557756p+0, 0x0p+0, 0x0p+0, 0x1.5a996b76709fap+1}, {{Gauss {0x1.9a79fec99f1aep-5, 0x1.8033daf8df7a5p+2, 0x1.a46822ff08894p+0}, Gauss {0x1.9f9acffa7eb6cp-5, 0x1.8030ced4e4c94p+2, 0x1.c8f42fe82517ep-1}}}}},
      {53, {Pack {Spd {-0x1.81d0ef1348b22p+6, -0x1.e8bb8f57f737ep+5, 0x0p+0}, Spd {-0x1.cfd0c3d25247dp+3, -0x1.7942d05f28848p+2, 0x0p+0}, Spd {0x1.c01098d477bcp+2, 0x1.3a2845996744bp+1, 0x0p+0}, Spd {}, 0x1.b438e086bdf4cp+3, 0x1.d273ffac1d29ep+2, 0x1.dfb167ec7863cp+3, 0x1.7dd99cbee807cp+2, 0x1.50a4fca42aed1p+1, 0x1.ff03978f28518p-1, 0x0p+0, 0x0p+0, 0x1.fd7cc39ffd60fp+0}, {{Gauss {-0x1.0d45e9185cee1p-3, 0x1.4d35efa615a8ep+2, 0x1.bfb2edfe75bc4p+0}, Gauss {-0x1.2e429e0a41a26p-5, 0x1.80a5c1c6088d7p+2, 0x1.5aed80a17b0f7p+1}}}}},
      {54, {Pack {Spd {0x1.7aeaedb4a3e12p+2, -0x1.5bee6b6177ea2p+6, 0x0p+0}, Spd {-0x1.cd6b8b69552e3p+1, -0x1.3de973cc8076bp+2, 0x0p+0}, Spd {0x1.3f5d74babcffdp+2, 0x1.58b1c864883fdp+1, 0x0p+0}, Spd {}, 0x1.17ff327aa68f5p+1, 0x1.8fce2e2b8c75cp+3, 0x1.39d8f92af9a8dp+2, 0x1.de7d41df9e60bp+3, 0x1.3d09833680c59p+1, 0x1.8e117095bf448p+2, 0x0p+0, 0x0p+0, 0x1.cb7ba6698bb4dp+0}, {{Gauss {-0x1.f6daf3003d3cap-3, 0x1.076cd0e3b2c26p+1, 0x1.b7a2f207efb89p+0}}}}},
      {7, {Pack {Spd {-0x1.8aaf74cd3176ap+5, -0x1.7c13f077ccc03p+5, 0x0p+0}, Spd {-0x1.c2002c0a4a05ep+3, -0x1.40b399f5dfeb9p+4, 0x0p+0}, Spd {0x1.0398958d9b5e9p+1, 0x1.28283d35eb745p+1, 0x0p+0}, Spd {}, 0x1.7cf403dddb121p+3, 0x1.782645e4e69fp+3, 0x1.d64ee392e1ef7p+2, 0x1.59d536933a04p+3, 0x1.22ff9f87f023fp+0, 0x1.2493804ff184fp+0, 0x0p+0, 0x0p+0, 0x1.6a4f4c6e6d9bep+1}, {{Gauss {0x1.806db50f40e5ap+0, 0x1.79ac68a936c59p+2, 0x1.b5f30e7ff583ap+0}, Gauss {-0x1.817a46173b85fp+0, 0x1.804c51116a8b9p+2, 0x1.b7558a7610279p+0}}}}},
      {56, {Pack {Spd {-0x1.434a3868265ap+3, -0x1.ab291aab7cf0fp+2, 0x0p+0}, Spd {-0x1.400009a59b2fap+3, -0x1.40001599c5388p+3, 0x0p+0}, Spd {0x1.ed02a9fe6bab5p+0, 0x1.73b5b1fe146d7p+0, 0x0p+0}, Spd {}, 0x1.35955a0458545p+2, 0x1.0fe9f77ffebdep+1, 0x1.98dc486ad2dcbp+1, 0x1.1b9207d4e0978p+1, 0x1.9985babf840efp-2, 0x1.6806734c27f6cp+1, 0x0p+0, 0x0p+0, 0x1.ffff1bd471dccp-2}, {{}}}},
      {9, {Pack {Spd {-0x1.b9bdc011d3672p+6, -0x1.a6bd7cf5f4e44p+6, 0x0p+0}, Spd {-0x1.833f5cf2495e1p+5, -0x1.bbea209aaa3adp+4, 0x0p+0}, Spd {0x1.2d58f7121ab4bp+2, 0x1.3edeebb341e15p+1, 0x0p+0}, Spd {}, 0x1.4fe4b2314013fp+3, 0x1.da26f60e0eb68p+3, 0x1.012dd4845133p+4, 0x1.cd6379b77c02bp+3, 0x1.749d5a187a4a5p-1, 0x1.4bd3369e4cb1bp+0, 0x0p+0, 0x0p+0, 0x1.adf11f926c7ebp+1}, {{Gauss {-0x1.8ea7ce0fc2ebap-7, 0x1.81823c85c24c4p+2, 0x1.db5b1c864884p+0}, Gauss {-0x1.75d13d74d594fp-9, 0x1.803ce63a5c1c6p+2, 0x1.516da0168b5ccp+1}}}}},
      {8, {Pack {Spd {-0x1.5bf8d5842b735p+6, -0x1.1f84b09e98dcep+6, 0x0p+0}, Spd {-0x1.699f077ccc038p+5, -0x1.8c0a4d2b2bfdbp+4, 0x0p+0}, Spd {0x1.e5f5275ee99a6p+1, 0x1.31d7ecbb7f9d7p+1, 0x0p+0}, Spd {}, 0x1.f82f2f9873ffbp+3, 0x1.b4edb2f661f19p+3, 0x1.53e08aefb2aaep+3, 0x1.8cfebaf102364p+3, 0x1.30116ebd4cfd1p-1, 0x1.ba219aaf6e422p-1, 0x0p+0, 0x0p+0, 0x1.9bc9ff92f2b67p+1}, {{Gauss {-0x1.21919ac79702ep+0, 0x1.8028954a7f801p+2, 0x1.9b78bbd380453p+0}, Gauss {0x1.234cd31769a91p+0, 0x1.7cd530489d27cp+2, 0x1.99306a2b1705p+0}}}}},
      {55, {Pack {Spd {-0x1.9a1169b4cf8p+1, -0x1.bec53b0813cacp+0, 0x0p+0}, Spd {-0x1.34ac2d1a147e4p-1, -0x1.7c122be70e977p+2, 0x0p+0}, Spd {0x1.cc4ab45938808p+1, 0x1.d9dd5687cc11cp-1, 0x0p+0}, Spd {}, 0x1.148d6b8ded93fp+1, 0x1.5a7fc4dc3d832p+2, 0x1.0aac075a6754ap+2, 0x1.929773ba0bfffp+2, 0x1.99253bf9ab523p-2, 0x1.9306aa3e2c8b8p+2, 0x0p+0, 0x0p+0, 0x1.0c33eae917863p-1}, {{}}}},
      {51, {Pack {Spd {-0x1.c375232d2bb23p+5, -0x1.d6f59253543afp+4, 0x0p+0}, Spd {-0x1.d96a39c51dabep+3, -0x1.68b28522ea0fdp+1, 0x0p+0}, Spd {0x1.2be8b3b320536p+1, 0x1.e665e02ea960bp+0, 0x0p+0}, Spd {}, 0x1.279ff7164c72ap+3, 0x1.9666666666666p+2, 0x1.51c58255b035cp+2, 0x1.9p+2, 0x1.3654d61b2a27fp+1, 0x1.790647c2b74a3p+0, 0x0p+0, 0x0p+0, 0x1.0463f9a49c2c2p+1}, {{Gauss {0x1.80427418d690ap+1, 0x1.805785f8d2e51p+2, 0x1.b4c447c30d307p-1}, Gauss {-0x1.3586ca89fc6dap-6, 0x1.80bc0e38a7e74p+2, 0x1.658b3700474dap+1}}}}},
      {4, {Pack {Spd {-0x1.143c6c97d8cf4p+4, -0x1.69bc5bd0e12e8p+3, 0x0p+0}, Spd {-0x1.fb248d7e02646p+1, -0x1.63ed740c4156ep+1, 0x0p+0}, Spd {0x1.c13faf42784a9p-1, 0x1.823dc486ad2ddp+0, 0x0p+0}, Spd {}, 0x1.2069468017119p+3, 0x1.83a8deb0fadf3p+2, 0x1.a4e071c53f39dp+2, 0x1.202ac1094a2bap+3, 0x1.16e02a77a2cedp-1, 0x1.82745bf26f1dcp+0, 0x0p+0, 0x0p+0, 0x1.97f1f9acffa7fp+0}, {{Gauss {0x1.a1aeb3dd11be7p+0, 0x1.56239e6ab9b24p+1, 0x1.caabef06b3786p+0}, Gauss {-0x1.0e33e78e1932dp+1, 0x1.f7f5c6c11a112p+0, 0x1.c180c308feac4p+0}}}}},
      {19, {Pack {Spd {-0x1.109e108c3f3ep+2, -0x1.523db99ef29efp+1, 0x0p+0}, Spd {-0x1.b33d519a26875p-2, -0x1.99937ecd6579p+1, 0x0p+0}, Spd {0x1.9ece6e8d7c54ep-1, 0x1.ea693e87fb0bap-1, 0x0p+0}, Spd {}, 0x1.b1d85e65e6e4ap+2, 0x1.bf8871532972cp+1, 0x1.2b1cea1b922ccp+3, 0x1.5eef5a964e8b8p+1, 0x1.a8c4552f067ddp+0, 0x1.00e82b2217acfp+1, 0x0p+0, 0x0p+0, 0x1.734f7c3f15bf8p-1}, {{}}}},
      {18, {Pack {Spd {0x1.11a9d769dec0fp+2, -0x1.2c66d1fbe0b69p+6, 0x0p+0}, Spd {-0x1.374bec679cc75p+1, -0x1.6fc566c1d5f8cp+4, 0x0p+0}, Spd {0x1.f6def2693e88p-1, 0x1.7ffb549f94856p+2, 0x0p+0}, Spd {}, 0x1.9e6e4f68f3df6p+1, 0x1.3ff9b604336b6p+4, 0x1.05b7e5cfd3119p+4, 0x1.5de5fb069bfb7p+3, 0x1.a242b22c37967p-1, 0x1.0cf11f3c86432p+2, 0x0p+0, 0x0p+0, 0x1.1b0cf0a94ae78p+1}, {{Gauss {-0x1.79ad4cd4c4e8dp-1, 0x1.f1b251208404p+1, 0x1.6e1d323fee2cap-1}}}}},
      {17, {Pack {Spd {-0x1.9281c9f72f76ep+6, -0x1.acea487336588p+5, 0x0p+0}, Spd {-0x1.b874fb549f948p+4, -0x1.730168b5cbff4p+3, 0x0p+0}, Spd {0x1.1f83cf2cf95d5p+1, 0x1.13544bb1af3a1p+1, 0x0p+0}, Spd {}, 0x1.0037b5aea3162p+4, 0x1.e16bf8769ec2dp+2, 0x1.018a2877ee4e2p+3, 0x1.e0440f238972p+2, 0x1.bd966be7afa72p+1, 0x1.b303292a45cecp-1, 0x0p+0, 0x0p+0, 0x1.4236c15d2d01cp+1}, {{Gauss {-0x1.5f6b1a2a4db16p-3, 0x1.800d23d4f15e8p+2, 0x1.166687f455a7dp+0}, Gauss {-0x1.b8fde2ef4e011p-7, 0x1.f77446f9b994ep+0, 0x1.257d73c925786p+1}}}}},
      {16, {Pack {Spd {-0x1.8f29b845564b6p+5, -0x1.6324028e4fb98p+5, 0x0p+0}, Spd {-0x1.1a7a97e132b56p+3, -0x1.02ecdf266ba49p+3, 0x0p+0}, Spd {0x1.e424b33daf8dfp+0, 0x1.a8b26394face6p+0, 0x0p+0}, Spd {}, 0x1.1ede8d5410f95p+3, 0x1.3efb3311a543fp+3, 0x1.b24cc6822ff09p+2, 0x1.fe1886df82b1fp+2, 0x1.02ad70e6f2e8cp+2, 0x1.84881bc2d4258p+0, 0x0p+0, 0x0p+0, 0x1.2285b9e8c47a1p+1}, {{Gauss {-0x1.98c586876e1dfp-2, 0x1.800af5fd47beep+2, 0x1.ec9b62c77574fp-1}, Gauss {-0x1.c1bb8c32a8c9cp-5, 0x1.801e3a7daa4fdp+2, 0x1.947735c182ecbp+0}}}}},
      {15, {Pack {Spd {-0x1.434e054690de1p+5, -0x1.d97d24180d3dp+4, 0x0p+0}, Spd {-0x1.93b547e06961cp+3, -0x1.0a3e186983516p+2, 0x0p+0}, Spd {0x1.023f811f4f50ap+1, 0x1.81361dc93ea2dp+0, 0x0p+0}, Spd {}, 0x1.f34da9003eea2p+2, 0x1.a79524bfd2e94p+2, 0x1.4bf6f8f041462p+2, 0x1.83f7d73c92578p+2, 0x1.8af587d6f9768p+0, 0x1.be740e19c929dp+0, 0x0p+0, 0x0p+0, 0x1.f0c6d612c6ac2p+0}, {{Gauss {-0x1.390c2c5e2cdcp-1, 0x1.ff4d37c1376d5p+0, 0x1.96d8f4f93bc0ap-1}, Gauss {-0x1.80c1fc8f32379p-4, 0x1.ff94855da2728p+0, 0x1.e92220bc382a1p+0}}}}},
      {38, {Pack {Spd {-0x1.5d6312bf1f652p+3, -0x1.fef2078458bfbp+2, 0x0p+0}, Spd {-0x1.40000ee3c891dp+3, -0x1.6b3c08c4ca9fdp+2, 0x0p+0}, Spd {0x1.478a3eb2c3365p+0, 0x1.6428f5c28f5c3p+0, 0x0p+0}, Spd {}, 0x1.42503968d759fp+2, 0x1.939fca3542d7fp+1, 0x1.f6d5d8e50903bp+1, 0x1.9f741eb5cdaccp+1, 0x1.84ed11e6b2512p-1, 0x1.59ce4735c9908p+1, 0x0p+0, 0x0p+0, 0x1.a591c30b82ebp+0}, {{}}}},
      {14, {Pack {Spd {-0x1.ac3739f340d4ep+4, -0x1.6d04a6223e187p+4, 0x0p+0}, Spd {-0x1.6e5ac471b4784p+1, -0x1.f77164c729f5ap+1, 0x0p+0}, Spd {0x1.a294467381d7ep+0, 0x1.50268900c521ep+0, 0x0p+0}, Spd {}, 0x1.4305425f2021p+2, 0x1.b099780baa583p+2, 0x1.7cbd5992428d4p+2, 0x1.4a52b0a6fc58bp+2, 0x1.d6f4384ba0e84p-1, 0x1.590c791fa68e7p+1, 0x0p+0, 0x0p+0, 0x1.11623076c050cp+1}, {{Gauss {-0x1.8ff972474538fp-2, 0x1.8000e27e0ef9ap+2, 0x1.43b7d84901d19p-1}, Gauss {0x1.d510d38cda6e7p-5, 0x1.8075afaf85943p+2, 0x1.028eef1bac2dfp+1}}}}},
      {37, {Pack {Spd {-0x1.25e4360a7e084p+2, -0x1.8186a16f50f0bp+1, 0x0p+0}, Spd {-0x1.82dce7cd03537p+3, -0x1.fffb33e864d8cp+0, 0x0p+0}, Spd {0x1.0000ae1049236p+2, 0x1.03720c8cd63ccp+0, 0x0p+0}, Spd {}, 0x1.28d2b579c9f08p+3, 0x1.abd23c785652fp+3, 0x1.400003e136106p+4, 0x1.300003ab862b2p+4, 0x1.3feaf0a4426ap+2, 0x1.7780a7b032b84p+0, 0x0p+0, 0x0p+0, 0x1.ff3fa383fe26ep-1}, {{}}}},
      {33, {Pack {Spd {-0x1.340f345069a4ep+5, -0x1.1938255b035bdp+5, 0x0p+0}, Spd {-0x1.076de54b48d3bp+3, -0x1.411cda2b5a20ep+2, 0x0p+0}, Spd {0x1.516e3f78bbd38p+1, 0x1.b43211cb039efp+0, 0x0p+0}, Spd {}, 0x1.193f7f06705c9p+3, 0x1.093126e978d5p+3, 0x1.59788db0574b4p+2, 0x1.06bb2788db057p+3, 0x1.f376f6d76252p+0, 0x1.8c4c19b4062c7p+0, 0x0p+0, 0x0p+0, 0x1.cb62d83c6c97ep+0}, {{Gauss {-0x1.d72324c836651p-2, 0x1.fbad6cb535009p+0, 0x1.163810e8858ffp+0}, Gauss {-0x1.6c8711d798d8bp-4, 0x1.fe3193f6c269ap+0, 0x1.11ed6ba8c5868p+1}}}}},
      {80, {Pack {Spd {-0x1.1c32170931013p+4, -0x1.254ac18f81e8ap+4, 0x0p+0}, Spd {-0x1.8cf9873ffac1dp+1, -0x1.bb655e28aa433p+1, 0x0p+0}, Spd {0x1.7a1522a6f3f53p+0, 0x1.3d6f08cc575cp+1, 0x0p+0}, Spd {}, 0x1.a7fb69984a0e4p+2, 0x1.d6b27243137bp+3, 0x1.54751efb6dcap+3, 0x1.000307f23cc8ep+4, 0x1.04a5d6bebe165p+1, 0x1.06e222b4e0c7dp+1, 0x0p+0, 0x0p+0, 0x1.87854046412cfp+0}, {{Gauss {0x1.152d234eb9a17p+0, 0x1.9fc842fa50939p+2, 0x1.31f51697f1f9bp+0}, Gauss {-0x1.8b7b28954a7f8p-4, 0x1.f6906034f3fd9p+1, 0x1.5046c764adff8p+1}}}}}
    }, {}
  };
  // clang-format on
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
