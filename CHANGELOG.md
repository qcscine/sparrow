# Changelog

## Release 2.0.0

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

## Release 1.0.1

Hotfix to allow compilation on OSX using Clang.

## Release 1.0.0

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
