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


#Calculate the dipole (IN DEBYE)
dip_test = scf.hf.dip_moment(ethanolPopleBasis, dm)
(sum(dip_test*dip_test))**(1/2)
