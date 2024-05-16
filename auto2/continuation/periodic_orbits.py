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

        # plots default behaviours
        self._default_marker = ''
        self._default_markersize = 6.
        self._default_linestyle = '-'
        self._default_linewidth = 1.2

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
                    return self.continuation[1].data[0].diagnostics[idx]['Multipliers']
                else:
                    warnings.warn('No backward branch to show the diagnostic for.')
                    return None
            else:
                idx = self.find_solution_index(idx)
                if idx is not None:
                    return self.continuation[0].data[0].diagnostics[idx]['Multipliers']
                else:
                    warnings.warn('No backward branch to show the diagnostic for.')
                    return None

        if idx >= 0:
            return self.continuation[0].data[0].diagnostics[idx]['Multipliers']
        else:
            if self.continuation[1] is not None:
                return self.continuation[1].data[0].diagnostics[-idx]['Multipliers']
            else:
                warnings.warn('No backward branch to show the diagnostic for.')
                return None

