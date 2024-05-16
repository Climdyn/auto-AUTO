import os
import sys
import warnings

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TABLEAU_COLORS

try:
    auto_directory = os.environ['AUTO_DIR']

    for path in sys.path:
        if auto_directory in path:
            break
    else:
        # sys.path.append(auto_directory + '/python/auto')
        sys.path.append(auto_directory + '/python')
except KeyError:
    warnings.warn('Unable to find auto directory environment variable.')

from auto2.parsers.config import ConfigParser
from auto2.continuation.fixed_points import FixedPointContinuation
from auto2.continuation.periodic_orbits import PeriodicOrbitContinuation


# TODO: add logging information

class BifurcationDiagram(object):

    def __init__(self, model_name):

        self.initial_points = None
        self.model_name = model_name
        self.config_object = ConfigParser('c.'+model_name)
        self.points_list = None
        self.solutions_list = None
        self.variables_list = None
        self.branches = dict()
        self.fp_computed = False
        self.po_computed = False

    def compute_fixed_points_diagram(self, initial_points=None, extra_comparison_parameters=None, **continuation_kwargs):

        if self.fp_computed:
            warnings.warn('Fixed point bifurcation diagram already computed. Aborting.')
            return None

        if initial_points is not None:
            self.initial_points = initial_points

        br_num = 1
        ncomp = 0
        for point in initial_points:
            parameters = point['parameters']
            initial_data = point['initial_data']

            used_continuation_kwargs = continuation_kwargs.copy()
            if 'PAR' not in used_continuation_kwargs:
                used_continuation_kwargs['PAR'] = dict()
            for par, val in parameters.items():
                if par not in used_continuation_kwargs['PAR']:
                    used_continuation_kwargs['PAR'][par] = val

            used_continuation_kwargs['IBR'] = br_num

            fp = FixedPointContinuation(model_name=self.model_name, config_object=self.config_object)
            fp.make_continuation(initial_data, **used_continuation_kwargs)

            for n, psol in self.branches.items():
                cpar = continuation_kwargs['ICP'][0]
                if isinstance(cpar, int):
                    cpar = self.config_object.parnames[cpar]
                if extra_comparison_parameters is not None:
                    cpar_list = extra_comparison_parameters.copy()
                    cpar_list.append(par)
                else:
                    cpar_list = [cpar]
                if fp.same_solutions_as(psol['continuation'], cpar_list):
                    warnings.warn('Not saving results of initial point '+str(ncomp)+' because it already exists (branch '+str(n)+').'
                                  '\nSkipping to next one.')  # should be a log instead
                    break
            else:
                self.branches[fp.branch_number] = {'parameters': parameters, 'continuation': fp, 'continuation_kwargs': used_continuation_kwargs}
                br_num += 1
            ncomp += 1

        self.fp_computed = True

    def plot_fixed_points_diagram(self, variables=(0, 1), ax=None, figsize=(10, 8), cmap=None, **kwargs):

        if 'plot_kwargs' not in kwargs:
            kwargs['plot_kwargs'] = dict()

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.gca()

        if cmap is not None:
            cmap = plt.get_cmap(cmap)

        handles = list()
        for i, b in enumerate(self.branches):
            kwargs['plot_kwargs']['label'] = 'BR '+str(b)
            if cmap is None:
                kwargs['plot_kwargs']['color'] = list(TABLEAU_COLORS.keys())[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_branches)
            self.branches[b]['continuation'].plot_branche_parts(variables, ax=ax, **kwargs)
            kwargs['plot_kwargs']['ls'] = '-'
            handles.append(Line2D([], [], **(kwargs['plot_kwargs'])))

        ax.legend(handles=handles)

    def plot_fixed_points_diagram_3D(self, variables=(0, 1), ax=None, figsize=(10, 8), cmap=None, **kwargs):

        if 'plot_kwargs' not in kwargs:
            kwargs['plot_kwargs'] = dict()

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(projection='3d')

        if cmap is not None:
            cmap = plt.get_cmap(cmap)

        handles = list()
        for i, b in enumerate(self.branches):
            kwargs['plot_kwargs']['label'] = 'BR '+str(b)
            if cmap is None:
                kwargs['plot_kwargs']['color'] = list(TABLEAU_COLORS.keys())[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_branches)
            self.branches[b]['continuation'].plot_branche_parts_3D(variables, ax=ax, **kwargs)
            kwargs['plot_kwargs']['ls'] = '-'
            handles.append(Line2D([], [], **(kwargs['plot_kwargs'])))

        ax.legend(handles=handles)
    @property
    def number_of_branches(self):
        return len(self.branches.keys())
