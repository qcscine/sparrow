Changelog
=========

Release 5.1.0
-------------

- Make compatible with Unity builds

Release 5.0.0
-------------

- Improve support for compilation on Windows (MSVC)
- Update address in license

Release 4.0.0
-------------

- Various bugfixes and improvements

Release 3.1.0
-------------

- Added parameter sets trans3d-0-1 (Sc, Ti, Fe, Co, Ni) and borg-0-1 (B)
- Enforce "C" locale for Molden files
- Various bugfixes and improvements

Release 3.0.1
-------------

- Various bugfixes and improvements

Release 3.0.0
-------------

- Excited-state calculations with NDDO-CIS and TD-DFTB
- Added a binary for the calculation of IR and UV/Vis spectra of trajectories
- Orbital steering calculations possible with binary (using a single system, single
  point) and Python bindings (can handle an individual system and systems along a
  trajectory)
- Can now set symmetry number for thermodynamics calculation
- Added density matrix RMSD SCF convergence check
- Calculate atomic second derivatives
- Removed Sparrow-specific Python bindings in favour of the more general ``Core::Calculator`` Python bindings
- Added Python bindings for excited-state calculators; see Utilities for more 
  infos on them
- Added patching functionality in DFTB embedded parameters, i.e. if the znorg-0-1 
  parameter set is chosen, then parameters are automatically sorted out between
  znorg-0-1 and mio-0-1 (patch parameters available: znorg-0-1, hyb-0-2)
- Fixed bug causing an instability in the calculation of the gradients (and Hessian matrix) in DFTBx
- Added ``conanfile.py`` for easier compilation and dependency-handling by Conan
- Made Sparrow relocatable by embedding the parameters in the compiled program;
  it is still possible to give parameter files externally
- Multiple bug fixes and stability fixes
- Improved testing

Release 2.0.1
-------------

- Corrected bug in setting the element collection in Python bindings
- Added missing include statement that prevented compilation on GCC 10.0.1

Release 2.0.0
-------------

- Calculate thermochemical properties
- Explicitly symmetrized Hessians
- Generation of molden input files for molecular orbital visualization
- General stability/performance improvements
- Access to the calculations and molden file generations also through SCINE Core::Calculator Python bindings
- Addition of 3ob-3-1 parameters sets for DFTB. Spin constants are extracted from
  Christof Köhler, Berücksichtigung von Spinpolarisationseffekten in einem dichtefunktionalbasierten Ansatz,
  PhD thesis, Departement Physik der Fakultät für Naturwissenschaften an der Universität Paderborn, 2004
- Addition of automatic linking to MKL/LAPACK/BLAS through Eigen
- Various bugfixes and improvements

Release 1.0.1
-------------

- Hotfix to allow compilation on OSX using Clang

Release 1.0.0
-------------

Initial release with the following features:

- Calculate electronic energies, nuclear gradients, and Hessians
- Calculate bond orders
- Restricted and unrestricted formalisms are implemented
- The following methods are available:
  - MNDO
  - AM1
  - RM1
  - PM3
  - PM6
  - non-SCC DFTB (DFTB0)
  - DFTB2
  - DFTB3
- All functionality can be accessed through a standalone binary
- Parallelized with OpenMP
