__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import setuptools

# Read README.rst for the long description
with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()


class EmptyListWithLength(list):
    """ Makes the wheel a binary distribution and platlib compliant. """

    def __len__(self):
        return 1


# Define the setup
setuptools.setup(
    name="scine_sparrow",
    version="@Sparrow_VERSION@",
    author="ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group",
    author_email="scine@phys.chem.ethz.ch",
    description="Open source semi-empirical quantum chemistry code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.scine.ethz.ch",
    packages=["scine_sparrow"],
    package_dir={"scine_sparrow": "scine_sparrow"},
    package_data={"scine_sparrow": ["*.module.so", "*.txt" @sparrow_PY_DEPS@]},
    install_requires=["scine_utilities"],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    zip_safe=False,
    test_suite='pytest',
    tests_require=['pytest'],
    ext_modules=EmptyListWithLength()
)
