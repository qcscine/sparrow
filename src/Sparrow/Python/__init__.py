__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import scine_utilities as utils
from distutils import ccompiler

shlib_suffix = ccompiler.new_compiler().shared_lib_extension
sparrow_module_filename = "sparrow.module" + shlib_suffix

manager = utils.core.ModuleManager.get_instance()
current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
if not manager.module_loaded('Sparrow'):
    if os.path.exists(os.path.join(current_path, sparrow_module_filename)):
        manager.load(os.path.join(current_path, sparrow_module_filename))
        if not manager.module_loaded('Sparrow'):
            raise ImportError("The {} was found but could not be loaded.".format(
                sparrow_module_filename))
    elif os.path.exists(os.path.join(lib_path, sparrow_module_filename)):
        manager.load(os.path.join(lib_path, sparrow_module_filename))
        if not manager.module_loaded('Sparrow'):
            raise ImportError("The {} was found but could not be loaded.".format(
                sparrow_module_filename))
    else:
        raise ImportError("The {} could not be located.".format(
            sparrow_module_filename))
