import os
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
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
        self.points_list = None  # TODO: to be implemented as property
        self.solutions_list = None  # TODO: to be implemented as property
        self.variables_list = None  # TODO: to be implemented as property
        self.branches = dict()
        self.fp_computed = False
        self.fp_parent = dict()
        self.fp_bp_computed = list()
        self.po_computed = False
        self._comparison_solutions_types = ('HB', 'BP', 'UZ', 'PD', 'EP', 'TR')

    def compute_fixed_points_diagram(self, initial_points=None, extra_comparison_parameters=None, comparison_tol=2.e-2, **continuation_kwargs):

        if self.fp_computed:
            warnings.warn('Fixed point bifurcation diagram already computed. Aborting.')
            return None

        if initial_points is not None:
            self.initial_points = initial_points

        if 'MNX' not in continuation_kwargs:
            continuation_kwargs['NMX'] = 9000
            warnings.warn('NMX parameters was not set, so setting it to 9000 points.')

        br_num = 1
        ncomp = 1
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

            valid_branch = self._check_continuation_against_other_branch(ncomp, fp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

            if valid_branch:
                self.branches[fp.branch_number] = {'parameters': parameters, 'continuation': fp, 'continuation_kwargs': used_continuation_kwargs}
                self.fp_parent[fp.branch_number] = None
                br_num += 1
            ncomp += 1

        bp_list = list()

        while True:
            nrecomp = 0

            new_branches = dict()

            for parent_branch_number, branch in self.branches.items():
                branching_points = branch['continuation'].get_filtered_solutions_list(labels='BP')

                if parent_branch_number not in self.fp_bp_computed:
                    used_continuation_kwargs = continuation_kwargs.copy()
                    used_continuation_kwargs['ISW'] = -1

                    for bp in branching_points:
                        for bpt in bp_list:
                            if self._check_if_solutions_are_close(bp, bpt, extra_comparison_parameters, comparison_tol):
                                break
                        else:
                            used_continuation_kwargs['IBR'] = br_num
                            fp = FixedPointContinuation(model_name=self.model_name, config_object=self.config_object)
                            fp.make_continuation(bp, **used_continuation_kwargs)

                            valid_branch = self._check_continuation_against_other_branch(ncomp, fp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                            if valid_branch:
                                new_branches[fp.branch_number] = {'parameters': bp.PAR, 'continuation': fp, 'continuation_kwargs': used_continuation_kwargs}
                                self.fp_parent[fp.branch_number] = parent_branch_number
                                br_num += 1
                            bp_list.append(bp)
                            ncomp += 1
                    self.fp_bp_computed.append(parent_branch_number)
                    nrecomp += 1

            self.branches.update(new_branches)

            if nrecomp == 0:
                break

        self.fp_computed = True

    @staticmethod
    def _check_if_solutions_are_close(sol1, sol2, extra_comparison_parameters, tol):

        cpar = sol1.c['ICP'][0]
        if extra_comparison_parameters is not None:
            comparison_parameters = [cpar]
            comparison_parameters.extend(extra_comparison_parameters)
        else:
            comparison_parameters = [cpar]

        npar = len(comparison_parameters)
        if isinstance(tol, float):
            tol = [tol] * npar

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        ssol1 = np.zeros(npar)
        ssol2 = np.zeros_like(ssol1)
        for i, params in enumerate(comparison_parameters):
            ssol1 = sol1[params]
            ssol2 = sol2[params]
        diff = ssol1 - ssol2
        return np.all(np.abs(diff) < tol)

    def _check_continuation_against_other_branch(self, ncomp, continuation, continuation_kwargs, extra_comparison_parameters, tol):

        fp = continuation
        initial_data = fp.initial_data

        valid_branch = False
        for n, psol in self.branches.items():
            cpar = continuation_kwargs['ICP'][0]
            if isinstance(cpar, int):
                cpar = self.config_object.parnames[cpar]
            if extra_comparison_parameters is not None:
                cpar_list = [cpar]
                cpar_list.extend(extra_comparison_parameters)
            else:
                cpar_list = [cpar]

            if fp.same_solutions_as(psol['continuation'], cpar_list, tol=tol, solutions_types=self._comparison_solutions_types):
                warnings.warn('Not saving results of initial point ' + str(ncomp) + ' because it already exists (branch ' + str(n) + ').'
                              '\nSkipping to next one.')  # should be a log instead
                break
            elif fp.solutions_in(psol['continuation'], cpar_list, tol=tol, solutions_types=self._comparison_solutions_types):
                warnings.warn('Not saving results of initial point ' + str(ncomp) + ' because it is already in branch ' + str(n) + '.'
                              '\nSkipping to next one.')  # should be a log instead
                break
            elif fp.solutions_part_of(psol['continuation'], cpar_list, tol=tol, forward=True, solutions_types=self._comparison_solutions_types):
                warnings.warn('Not storing full results of initial point ' + str(ncomp) + ' because it merges with branch ' + str(n) + '.'
                              '\nSaving only the relevant part.')  # should be a log instead
                _, common_solutions = fp.solutions_part_of(psol['continuation'], cpar_list, tol=tol,
                                                           return_solutions=True, forward=True, solutions_types=self._comparison_solutions_types)
                first_sol = common_solutions[0]
                continuation_kwargs['NMX'] = first_sol['PT'] + 1
                fp.make_forward_continuation(initial_data, **continuation_kwargs)
                valid_branch = True
            elif fp.solutions_part_of(psol['continuation'], cpar_list, tol=tol, forward=False):
                warnings.warn('Not storing full results of initial point ' + str(ncomp) + ' because it merges with branch ' + str(n) + '.'
                              '\nSaving only the relevant part.')  # should be a log instead
                _, common_solutions = fp.solutions_part_of(psol['continuation'], cpar_list, tol=tol,
                                                           return_solutions=True, forward=False, solutions_types=self._comparison_solutions_types)
                first_sol = common_solutions[0]
                continuation_kwargs['NMX'] = first_sol['PT'] + 1
                fp.make_backward_continuation(initial_data, **continuation_kwargs)
                valid_branch = True
            elif fp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol, forward=True):
                warnings.warn('Not storing full results of initial point ' + str(ncomp) + ' because it connects to branch ' + str(n) + '.'
                              '\nSaving only the relevant part.')  # should be a log instead
                _, sol = fp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol,
                                                  return_solutions=True, forward=True, solutions_types=self._comparison_solutions_types)
                continuation_kwargs['NMX'] = sol['PT'] + 1
                fp.make_forward_continuation(initial_data, **continuation_kwargs)
                valid_branch = True
                break
            elif fp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol, forward=False):
                warnings.warn('Not storing full results of initial point ' + str(ncomp) + ' because it connects to branch ' + str(n) + '.'
                              '\nSaving only the relevant part.')  # should be a log instead
                _, sol = fp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol,
                                                  return_solutions=True, forward=False, solutions_types=self._comparison_solutions_types)
                continuation_kwargs['NMX'] = sol['PT'] + 1
                fp.make_backward_continuation(initial_data, **continuation_kwargs)
                valid_branch = True
                break
        else:
            valid_branch = True

        return valid_branch

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
            kwargs['plot_kwargs']['linestyle'] = '-'
            handles.append(Line2D([], [], **(kwargs['plot_kwargs'])))

        ax.legend(handles=handles)

    def plot_fixed_points_diagram_3D(self, variables=(0, 1, 2), ax=None, figsize=(10, 8), cmap=None, **kwargs):

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
            kwargs['plot_kwargs']['linestyle'] = '-'
            handles.append(Line2D([], [], **(kwargs['plot_kwargs'])))

        ax.legend(handles=handles)

    @property
    def number_of_branches(self):
        return len(self.branches.keys())
