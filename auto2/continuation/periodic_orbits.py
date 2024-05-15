import os
import sys
import warnings

import matplotlib.pyplot as plt

auto_directory = os.environ['AUTO_DIR']
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

import auto.AUTOCommands as ac
import auto.runAUTO as ra
from auto.parseS import AUTOSolution
from auto2.continuation.base import Continuation


class PeriodicOrbitContinuation(Continuation):

    def __init__(self, model_name, config_object):

        Continuation.__init__(self, model_name, config_object)

    def make_continuation(self, initial_data, store_name="", only_forward=True, **continuation_kwargs):
        runner = ra.runAUTO()
        ac.load(self.model_name, runner=runner)

        if 'MXBF' in continuation_kwargs:
            warnings.warn('Disabling automatic continuation of branch points (MXBF set to 0)')
        continuation_kwargs['MXBF'] = 0

        self.initial_data = initial_data

        if isinstance(initial_data, AUTOSolution):
            cf = ac.run(initial_data, runner=runner, **continuation_kwargs)
            if not only_forward:
                if 'DS' in continuation_kwargs:
                    continuation_kwargs['DS'] = - continuation_kwargs['DS']
                    cb = ac.run(initial_data, runner=runner, **continuation_kwargs)
                else:
                    cb = ac.run(initial_data, DS='-', runner=runner, **continuation_kwargs)
            else:
                cb = None

        else:
            cf = ac.run(self.model_name, dat=initial_data, runner=runner, **continuation_kwargs)
            if not only_forward:
                if 'DS' in continuation_kwargs:
                    continuation_kwargs['DS'] = - continuation_kwargs['DS']
                    cb = ac.run(self.model_name, dat=initial_data, runner=runner, **continuation_kwargs)
                else:
                    cb = ac.run(self.model_name, DS='-', dat=initial_data, runner=runner, **continuation_kwargs)
            else:
                cb = None

        self.continuation = list([cf, cb])
        self.branch_number = self.continuation[0].data[0].BR

        if store_name:
            self.auto_save(store_name)

    def orbit_stability(self, idx):
        if isinstance(idx, str):
            if idx[0] == '-':
                idx = self.find_solution_index(idx)
                if idx is not None:
                    return self.continuation[1].data[0].diagnostics[idx]['Eigenvalues']
                else:
                    warnings.warn('No backward branch to show the diagnostic for.')
                    return None
            else:
                idx = self.find_solution_index(idx)
                if idx is not None:
                    return self.continuation[0].data[0].diagnostics[idx]['Eigenvalues']
                else:
                    warnings.warn('No backward branch to show the diagnostic for.')
                    return None

        if idx >= 0:
            return self.continuation[0].data[0].diagnostics[idx]
        else:
            if self.continuation[1] is not None:
                return self.continuation[1].data[0].diagnostics[-idx]
            else:
                warnings.warn('No backward branch to show the diagnostic for.')
                return None

    def plot_solutions(self, variables=(0, 1), ax=None, figsize=(10, 8), markersize=12., marker='', linestyle='-',
                       linewidth=1.2, plot_kwargs=None, labels=None, indices=None, parameter=None, value=None,
                       variables_name=None, tol=0.01):

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.gca()

        if plot_kwargs is None:
            plot_kwargs = dict()

        solutions_list = list()
        if self.continuation:
            solutions_list = self.solutions_list

            if labels is not None:
                if not isinstance(labels, (list, tuple)):
                    labels = [labels]
                new_solutions_list = list()
                for sol in solutions_list:
                    if sol['TY'] in labels:
                        new_solutions_list.append(sol)
                solutions_list = new_solutions_list
            elif indices is not None:
                if not isinstance(indices, (list, tuple)):
                    indices = [indices]
                new_solutions_list = list()
                for sol in solutions_list:
                    if sol['PT'] - 1 in indices:
                        new_solutions_list.append(sol)
                solutions_list = new_solutions_list
            elif parameter is not None:
                if value is not None:
                    new_solutions_list = list()
                    for sol in solutions_list:
                        if abs(sol[parameter] - value) < tol:
                            new_solutions_list.append(sol)
                    solutions_list = new_solutions_list

                else:
                    solutions_list = list()

        keys = self.config_object.variables

        if variables[0] in keys:
            var1 = variables[0]
        else:
            try:
                var1 = keys[variables[0]]
            except:
                var1 = keys[0]

        if variables[1] in keys:
            var2 = variables[1]
        else:
            try:
                var2 = keys[variables[1]]
            except:
                var2 = keys[1]

        for sol in solutions_list:
            x = sol[var1]
            y = sol[var2]
            line_list = ax.plot(x, y, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth,
                                **plot_kwargs)
            c = line_list[0].get_color()
            plot_kwargs['color'] = c

        if variables_name is None:
            ax.set_xlabel(var1)
            ax.set_ylabel(var2)
        else:
            if isinstance(variables_name, dict):
                ax.set_xlabel(variables_name[var1])
                ax.set_ylabel(variables_name[var2])
            else:
                ax.set_xlabel(variables_name[0])
                ax.set_ylabel(variables_name[1])
        return ax
