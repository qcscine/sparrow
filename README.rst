SCINE - Sparrow
===============

Introduction
------------

Sparrow is a code for fast semiempirical quantum chemical calculations. It
provides the following methods:

- MNDO
- AM1
- RM1
- PM3
- PM6
- non-SCC DFTB (DFTB0)
- DFTB2
- DFTB3

Sparrow can calculate electronic energies, nuclear gradients and
Hessians for the electronic ground state, as well as electronic vertical
transition energies and the electronic transition dipoles.

License and Copyright Information
---------------------------------

Sparrow is distributed under the BSD 3-clause "New" or "Revised" License.
For more license and copyright information, see the file ``LICENSE.txt`` in the
repository.

Installation and Usage
----------------------

For instructions on how to install and use Sparrow as well as for a detailed
documentation of the entire functionality of Sparrow, please consult the user
manual found in the ``manual`` directory in the repository.
Alternatively the manual can also be found on the official GitHub website,
SCINE website and in the hosted documentation.

How to Cite
-----------

When publishing results obtained with Sparrow, please cite the corresponding
release as archived on `Zenodo <https://zenodo.org/record/3244105>`_ (DOI
10.5281/zenodo.3244105; please use the DOI of the respective release).

In addition, we kindly request you to cite the following article when using Sparrow:

T. Husch, A. C. Vaucher, M. Reiher, "Semiempirical molecular orbital models
based on the neglect of diatomic differential overlap approximation", *Int.
J. Quantum Chem.*, **2018**, *118*, e25799.

Support and Contact
-------------------

In case you should encounter problems or bugs, please write a short message
to scine@phys.chem.ethz.ch.

Third-Party Libraries Used
--------------------------

SCINE Sparrow makes use of the following third-party libraries:

- `Boost <https://www.boost.org/>`_
- `Cereal <https://uscilab.github.io/cereal/>`_
- `Eigen <http://eigen.tuxfamily.org>`_
- `Google Test <https://github.com/google/googletest>`_
- `pybind11 <https://github.com/pybind/pybind11>`_
- `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_
