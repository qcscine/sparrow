__copyright__ = """This file is part of SCINE Sparrow.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from conans import ConanFile


class TestPackageConan(ConanFile):
    def build(self):
        pass

    def test(self):
        if self.options["scine_sparrow"].python:
            self.output.info("Trying to import 'scine_sparrow'")
            import scine_sparrow
            self.output.info("Import worked")
