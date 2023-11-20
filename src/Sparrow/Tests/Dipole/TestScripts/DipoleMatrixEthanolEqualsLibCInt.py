__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from pyscf import gto, scf
ethanolPopleBasis = gto.Mole()

ethanolPopleBasis.atom = '''
                  C       -0.025722743    -0.038148774     0.007736703;
                  C        1.494872990    -0.038148774     0.007736703;
                  O        2.017472251     1.299308464     0.007736703;
                  H       -0.418489864     0.516363236     0.871446070;
                  H       -0.435916359     0.434351771    -0.891900626;
                  H       -0.429207593    -1.056100478     0.060014382;
                  H        1.909273199    -0.583326143    -0.859094289;
                  H        1.916182172    -0.455131719     0.942989512;
                  H        1.715886732     1.789449608    -0.779284087
'''
ethanolPopleBasis.basis = 'sto-6g'
ethanolPopleBasis.build()

#Run the calculation
m2 = scf.hf.SCF(ethanolPopleBasis)
m2.kernel()
dm = m2.make_rdm1()
print(m2.mol.atom)


###
##
## PySCF source code was modified to output the Atomic Orbitals dipole operator integral matrix
## for the valence orbitals. Core orbitals and core charge was neglected
##
###

#AO dip in x:
#[[-4.86089394e-02 -1.18230573e-02  7.58402264e-02 -7.25567384e-20
#  -2.96164781e-20  4.82124045e-06  3.23851408e-03 -5.34339050e-03
#  -2.32713293e-22  3.33078177e-23  1.21176674e-14  8.07557567e-05
#  -1.07894766e-04 -8.54046121e-05 -5.35408171e-23 -6.35459337e-03
#  -6.55509492e-03 -6.49139432e-03  4.44040972e-04  4.38292103e-04
#   9.13729639e-05]
# [-1.18230573e-02 -4.86089394e-02  8.39170413e-01  0.00000000e+00
#  -4.10291795e-21  7.56525578e-02  4.30703630e-01 -2.89204265e-01
#   0.00000000e+00 -3.08034404e-20  5.49737108e-03  8.00293844e-02
#  -6.10088314e-02 -5.38500946e-02 -7.55966770e-20 -2.13523276e-01
#  -2.22815808e-01 -2.19448099e-01  1.63202848e-01  1.61991687e-01
#   6.18712270e-02]
# [ 7.58402264e-02  8.39170413e-01 -4.86089394e-02  0.00000000e+00
#   1.01928333e-19  1.29079567e-01  6.56204527e-01 -4.46310298e-01
#   0.00000000e+00  6.52969731e-20  7.89070318e-03  1.12091149e-01
#  -7.79785659e-02 -8.49717355e-02 -8.26469779e-20  4.20417154e-01
#   4.28236506e-01  4.25634018e-01  2.39870310e-01  2.38147404e-01
#   8.10457250e-02]
# [-7.25567384e-20  0.00000000e+00  0.00000000e+00 -4.86089394e-02
#   0.00000000e+00  8.72966528e-20  0.00000000e+00  0.00000000e+00
#   2.35428080e-01  0.00000000e+00  5.14077513e-03  6.32875150e-02
#  -5.76790972e-02 -2.33203514e-02 -4.80248628e-20 -1.00115253e-01
#  -8.91386839e-02  1.89110615e-01 -4.86094603e-02 -3.68282874e-02
#   5.36523023e-02]
# [-2.96164781e-20 -4.10291795e-21  1.01928333e-19  0.00000000e+00
#  -4.86089394e-02 -8.85359276e-21 -3.96788887e-20 -1.73435488e-20
#   0.00000000e+00  2.35428080e-01  3.08420145e-22 -3.88534017e-20
#   2.03553435e-20  1.93995391e-20  2.56527027e-02 -1.55939781e-01
#   1.69719355e-01 -9.71191856e-03 -7.72889505e-02  8.26023215e-02
#  -2.31043525e-02]
# [ 4.82124045e-06  7.56525578e-02  1.29079567e-01  8.72966528e-20
#  -8.85359276e-21  2.82490054e+00  6.87095036e-01  7.58402264e-02
#   4.21662295e-18  1.72115677e-18  3.14737509e-06  5.11428994e-02
#  -3.04052145e-02 -8.13727225e-02 -1.59249351e-20  1.42996384e-02
#   1.39474831e-02  1.40949074e-02  1.79701948e-01  1.78696204e-01
#   2.24281384e-02]
# [ 3.23851408e-03  4.30703630e-01  6.56204527e-01  0.00000000e+00
#  -3.96788887e-20  6.87095036e-01  2.82490054e+00  8.39170413e-01
#   0.00000000e+00  2.38440404e-19  8.31464739e-02  8.79306724e-01
#  -1.41410084e-01 -7.95526771e-01 -4.28394609e-19  1.06297773e-01
#   1.02987483e-01  1.04317566e-01  1.56319955e+00  1.56168420e+00
#   3.97001974e-01]
# [-5.34339050e-03 -2.89204265e-01 -4.46310298e-01  0.00000000e+00
#  -1.73435488e-20  7.58402264e-02  8.39170413e-01  2.82490054e+00
#   0.00000000e+00  1.01928333e-19  5.26852600e-02  5.52164041e-01
#   2.81883313e-01 -6.11679580e-01 -3.41229383e-20 -4.82149672e-02
#  -4.64788538e-02 -4.71871221e-02  9.06153855e-01  9.13091550e-01
#   1.44658129e-01]
# [-2.32713293e-22  0.00000000e+00  0.00000000e+00  2.35428080e-01
#   0.00000000e+00  4.21662295e-18  0.00000000e+00  0.00000000e+00
#   2.82490054e+00  0.00000000e+00  1.32447899e-01  1.04470546e+00
#  -2.82872611e-01 -7.82052816e-01 -1.72295830e-20  3.41021592e-02
#   2.80646489e-02 -6.13265728e-02 -7.35541922e-01 -5.61321730e-01
#   4.31570957e-01]
# [ 3.33078177e-23 -3.08034404e-20  6.52969731e-20  0.00000000e+00
#   2.35428080e-01  1.72115677e-18  2.38440404e-19  1.01928333e-19
#   0.00000000e+00  2.82490054e+00  2.27382821e-20 -9.54731645e-20
#  -9.86571881e-20 -1.82388254e-19  4.44199185e-01  5.31176130e-02
#  -5.34348713e-02  3.14947249e-03 -1.16951027e+00  1.25899088e+00
#  -1.85847897e-01]
# [ 1.21176674e-14  5.49737108e-03  7.89070318e-03  5.14077513e-03
#   3.08420145e-22  3.14737509e-06  8.31464739e-02  5.26852600e-02
#   1.32447899e-01  2.27382821e-20  3.81247002e+00  8.81232912e-01
#   5.37365790e-02 -3.65095974e-16  4.52807877e-19  3.61814097e-03
#   3.21215166e-03  7.07585503e-04  1.57291510e-02  1.91900131e-02
#   2.01557787e-01]
# [ 8.07557567e-05  8.00293844e-02  1.12091149e-01  6.32875150e-02
#  -3.88534017e-20  5.11428994e-02  8.79306724e-01  5.52164041e-01
#   1.04470546e+00 -9.54731645e-20  8.81232912e-01  3.81247002e+00
#   6.41499160e-01 -5.67844538e-18  4.95861631e-19  4.50780104e-02
#   4.04490875e-02  1.06234278e-02  2.56115903e-01  3.02908044e-01
#   1.65976414e+00]
# [-1.07894766e-04 -6.10088314e-02 -7.79785659e-02 -5.76790972e-02
#   2.03553435e-20 -3.04052145e-02 -1.41410084e-01  2.81883313e-01
#  -2.82872611e-01 -9.86571881e-20  5.37365790e-02  6.41499160e-01
#   3.81247002e+00 -5.64459277e-18  2.28697711e-19 -3.13013216e-02
#  -2.77180344e-02 -5.83894914e-03  3.15818830e-02  3.75785291e-02
#  -1.28533357e-01]
# [-8.54046121e-05 -5.38500946e-02 -8.49717355e-02 -2.33203514e-02
#   1.93995391e-20 -8.13727225e-02 -7.95526771e-01 -6.11679580e-01
#  -7.82052816e-01 -1.82388254e-19 -3.65095974e-16 -5.67844538e-18
#  -5.64459277e-18  3.81247002e+00  9.48612683e-34 -1.42829902e-02
#  -1.39328119e-02 -8.50582161e-03 -2.46156918e-01 -2.80408893e-01
#   6.92419791e-01]
# [-5.35408171e-23 -7.55966770e-20 -8.26469779e-20 -4.80248628e-20
#   2.56527027e-02 -1.59249351e-20 -4.28394609e-19 -3.41229383e-20
#  -1.72295830e-20  4.44199185e-01  4.52807877e-19  4.95861631e-19
#   2.28697711e-19  9.48612683e-34  3.81247002e+00  1.57563416e-02
#  -1.44914512e-02  1.88784463e-04 -1.13339277e-01  1.49479707e-01
#  -1.11182009e+00]
# [-6.35459337e-03 -2.13523276e-01  4.20417154e-01 -1.00115253e-01
#  -1.55939781e-01  1.42996384e-02  1.06297773e-01 -4.82149672e-02
#   3.41021592e-02  5.31176130e-02  3.61814097e-03  4.50780104e-02
#  -3.13013216e-02 -1.42829902e-02  1.57563416e-02 -7.90831229e-01
#  -1.39819431e-01 -1.37837229e-01  2.53977883e-02  7.02822613e-02
#   2.73933878e-02]
# [-6.55509492e-03 -2.22815808e-01  4.28236506e-01 -8.91386839e-02
#   1.69719355e-01  1.39474831e-02  1.02987483e-01 -4.64788538e-02
#   2.80646489e-02 -5.34348713e-02  3.21215166e-03  4.04490875e-02
#  -2.77180344e-02 -1.39328119e-02 -1.44914512e-02 -1.39819431e-01
#  -8.23762532e-01 -1.40886863e-01  6.60058244e-02  2.47223469e-02
#   5.84654050e-02]
# [-6.49139432e-03 -2.19448099e-01  4.25634018e-01  1.89110615e-01
#  -9.71191856e-03  1.40949074e-02  1.04317566e-01 -4.71871221e-02
#  -6.13265728e-02  3.14947249e-03  7.07585503e-04  1.06234278e-02
#  -5.83894914e-03 -8.50582161e-03  1.88784463e-04 -1.37837229e-01
#  -1.40886863e-01 -8.11084801e-01  6.63039131e-02  6.43047760e-02
#   7.69348038e-03]
# [ 4.44040972e-04  1.63202848e-01  2.39870310e-01 -4.86094603e-02
#  -7.72889505e-02  1.79701948e-01  1.56319955e+00  9.06153855e-01
#  -7.35541922e-01 -1.16951027e+00  1.57291510e-02  2.56115903e-01
#   3.15818830e-02 -2.46156918e-01 -1.13339277e-01  2.53977883e-02
#   6.60058244e-02  6.63039131e-02  3.60800344e+00  5.87482133e-01
#   2.18864769e-01]
# [ 4.38292103e-04  1.61991687e-01  2.38147404e-01 -3.68282874e-02
#   8.26023215e-02  1.78696204e-01  1.56168420e+00  9.13091550e-01
#  -5.61321730e-01  1.25899088e+00  1.91900131e-02  3.02908044e-01
#   3.75785291e-02 -2.80408893e-01  1.49479707e-01  7.02822613e-02
#   2.47223469e-02  6.43047760e-02  5.87482133e-01  3.62105951e+00
#   9.97649629e-02]
# [ 9.13729639e-05  6.18712270e-02  8.10457250e-02  5.36523023e-02
#  -2.31043525e-02  2.24281384e-02  3.97001974e-01  1.44658129e-01
#   4.31570957e-01 -1.85847897e-01  2.01557787e-01  1.65976414e+00
#  -1.28533357e-01  6.92419791e-01 -1.11182009e+00  2.73933878e-02
#   5.84654050e-02  7.69348038e-03  2.18864769e-01  9.97649629e-02
#   3.24255598e+00]]
#pyscf nuclear dipole = [55.20345014 19.99762741 -0.94693227]
#Dipole moment(X, Y, Z, Debye): -0.99717, -0.57778, -1.03090
