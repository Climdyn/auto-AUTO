
"""

Configuration object definition
===============================

This module defines the configuration object used in auto-AUTO.

"""

import os
import sys
import logging

# putting the logger creation here since this module is always called
logger = logging.getLogger('logger')
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s %(levelname)s: Module %(filename)s -- %(message)s')

fh = logging.FileHandler('auto2.log', mode='w', encoding='utf-8')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

try:
    auto_directory = os.environ['AUTO_DIR']

    for path in sys.path:
        if auto_directory in path:
            break
    else:
        # sys.path.append(auto_directory + '/python/auto')
        sys.path.append(auto_directory + '/python')
except KeyError:
    logger.warning('Unable to find auto directory environment variable.')

import auto.parseC as parseC


class ConfigParser(object):
    """An object to load and parse AUTO configuration files.

    Parameters
    ----------
    config_file: str
        Path to the |AUTO| config file.

    Attributes
    ----------
    config_path: str
        Path to the |AUTO| config file.
    config_object: AUTO Configuration object
        The native AUTO configuration Python object.

    """

    def __init__(self, config_file):

        self.config_path = config_file
        self.config_object = parseC.parseC(config_file)

    def keys(self):
        return self.config_object.keys()

    def __getitem__(self, item):
        return self.config_object.__getitem__(item)

    def __str__(self):
        return self.config_object.__str__()

    def __repr__(self):
        return self.config_object.__repr__()

    @property
    def variables(self):
        """list(str): List the names of the variables of the dynamical system."""
        variables_list = self['unames']
        return [v for n, v in variables_list]

    @property
    def parameters(self):
        """list(str): List the names of the available continuations parameters."""
        parameters_list = self['parnames']
        return [p for n, p in parameters_list]

    @property
    def variables_dict(self):
        """dict(str): Dictionary of the names of the variables of the dynamical system,
        indexed by their AUTO number."""
        variables_list = self['unames']
        return {n: v for n, v in variables_list}

    @property
    def parameters_dict(self):
        """dict(str): Dictionary of the names of the variables of the dynamical system,
        indexed by their AUTO number."""
        parameters_list = self['parnames']
        return {n: p for n, p in parameters_list}

    @property
    def parnames(self):
        """list: Link to AUTO parnames structure."""
        return self['parnames']

    @property
    def unames(self):
        """list: Link to AUTO unames structure."""
        return self['unames']

    @property
    def ndim(self):
        """int: Dimension of the dynamical system."""
        return self['NDIM']

    @property
    def continuation_parameters(self):
        """list(str): Parameters used by default for continuation."""
        return self['ICP']

    @property
    def parameters_solution_points(self):
        """list: List of user defined points parameters values."""
        return self['UZR']

    @property
    def parameters_bounds(self):
        """list: List of user defined parameters bound values."""
        return self['UZSTOP']

