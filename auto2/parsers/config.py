
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
    """
    An object to load and parse AUTO configuration files.
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
        variables_list = self['unames']
        return [v for n, v in variables_list]

    @property
    def parameters(self):
        parameters_list = self['parnames']
        return [p for n, p in parameters_list]

    @property
    def variables_dict(self):
        variables_list = self['unames']
        return {n: v for n, v in variables_list}

    @property
    def parameters_dict(self):
        parameters_list = self['parnames']
        return {n: p for n, p in parameters_list}

    @property
    def parnames(self):
        return self['parnames']

    @property
    def unames(self):
        return self['unames']

    @property
    def ndim(self):
        return self['NDIM']

    @property
    def continuation_parameters(self):
        return self['ICP']

    @property
    def parameters_solution_points(self):
        return self['UZR']

    @property
    def parameters_bounds(self):
        return self['UZSTOP']

