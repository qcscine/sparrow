import sys
from dev.conan.base import ScineConan


class ScineSparrowConan(ScineConan):
    name = "scine_sparrow"
    version = "5.1.0"
    url = "https://github.com/qcscine/sparrow"
    description = """Sparrow is a code for fast semiemprical quantum chemical calculations.
It provides the methods such as MNDO, AM1, RM1, PM3, PM6, DFTB0, DFTB2, and DFTB3.
Sparrow can calculate electronic energies as well as nuclear gradients and Hessians."""
    options = {
        "python": [True, False],
        "tests": [True, False],
        "coverage": [True, False],
        "microarch": ["detect", "none"],
        "python_version": "ANY"
    }
    python_version_string = str(sys.version_info.major) + \
                    "." + str(sys.version_info.minor)
    default_options = {
        "python": False,
        "tests": False,
        "coverage": False,
        "microarch": "none",
        "python_version": python_version_string
    }
    exports = "dev/conan/*.py"
    exports_sources = ["dev/cmake/*", "src/*", "CMakeLists.txt",
                       "README.rst", "LICENSE.txt", "dev/conan/hook.cmake",
                       "dev/conan/glue/*"]
    build_requires = "cereal/1.3.0"
    requires = "scine_utilities/[=10.0.0]"
    cmake_name = "Sparrow"

    def package_info(self):
        super().package_info()
        self.cpp_info.components["Sparrow"].libs = ["sparrow.module.so"]
