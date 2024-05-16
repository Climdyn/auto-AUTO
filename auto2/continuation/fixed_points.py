import os
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np

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


class FixedPointContinuation(Continuation):

    def __init__(self, model_name, config_object):

        Continuation.__init__(self, model_name, config_object)

        # plots default behaviours
        self._default_marker = 'x'
        self._default_markersize = 6.
        self._default_linestyle = ' '
        self._default_linewidth = 1.2

    def make_continuation(self, initial_data, store_name="", only_forward=False, **continuation_kwargs):
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
            u = {i + 1: initial_data[i] for i in range(self.config_object.ndim)}
            cf = ac.run(self.model_name, U=u, runner=runner, **continuation_kwargs)
            if not only_forward:
                if 'DS' in continuation_kwargs:
                    continuation_kwargs['DS'] = - continuation_kwargs['DS']
                    cb = ac.run(self.model_name, U=u, runner=runner, **continuation_kwargs)
                else:
                    cb = ac.run(self.model_name, DS='-', U=u, runner=runner, **continuation_kwargs)
            else:
                cb = None

        self.continuation = list([cf, cb])
        self.branch_number = self.continuation[0].data[0].BR

        if store_name:
            self.auto_save(store_name)

    def point_stability(self, idx):
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
            return self.continuation[0].data[0].diagnostics[idx]['Eigenvalues']
        else:
            if self.continuation[1] is not None:
                return self.continuation[1].data[0].diagnostics[-idx]['Eigenvalues']
            else:
                warnings.warn('No backward branch to show the diagnostic for.')
                return None

    def same_solutions_as(self, other, parameter, solutions_type=('HB', 'BP', 'UZ', 'PD'), tol=2.e-2):
        ssol = np.array(self.solutions_parameters(parameter, solutions_type))
        osol = np.array(other.solutions_parameters(parameter, solutions_type))

        osol.sort()
        ssol.sort()

        if len(ssol) != len(osol):
            return False
        else:
            dif = ssol - osol
            return np.all(np.abs(dif) < tol)



