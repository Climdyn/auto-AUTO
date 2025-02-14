import os
import sys
import warnings
import logging
import traceback
import glob

logger = logging.getLogger('logger')

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

import auto.AUTOCommands as ac
import auto.runAUTO as ra
from auto.parseS import AUTOSolution
from auto.AUTOExceptions import AUTORuntimeError
from auto2.continuations.base import Continuation

import auto2.continuations.periodic_orbits as poc


class FixedPointContinuation(Continuation):

    def __init__(self, model_name, config_object, path_name=None):

        Continuation.__init__(self, model_name, config_object, path_name)

        # plots default behaviours
        self._default_marker = 'x'
        self._default_markersize = 6.
        self._default_linestyle = ' '
        self._default_linewidth = 1.2

    def make_continuation(self, initial_data, auto_suffix="", only_forward=False, **continuation_kwargs):

        self.make_forward_continuation(initial_data, "", **continuation_kwargs)
        if not only_forward:
            self.make_backward_continuation(initial_data, "", **continuation_kwargs)

        if auto_suffix:
            self.auto_save(auto_suffix)

    def make_forward_continuation(self, initial_data, auto_suffix="", **continuation_kwargs):
        runner = ra.runAUTO()
        ac.load(self.model_name, runner=runner)

        if 'MXBF' in continuation_kwargs:
            warnings.warn('Disabling automatic continuation of branch points (MXBF set to 0)')
        continuation_kwargs['MXBF'] = 0

        if 'IBR' not in continuation_kwargs and self.branch_number is not None:
            continuation_kwargs['IBR'] = self.branch_number

        self.initial_data = initial_data

        if isinstance(initial_data, AUTOSolution):
            for retry in range(self._retry):
                try:
                    cf = ac.run(initial_data, runner=runner, **continuation_kwargs)
                except AUTORuntimeError:
                    print(traceback.format_exc())
                    warnings.warn('AUTO continuation failed, possibly retrying.')
                else:
                    break
            else:
                warnings.warn('Problem to complete the forward AUTO continuation, returning nothing.')
                cf = None

        else:
            u = {i + 1: initial_data[i] for i in range(self.config_object.ndim)}
            for retry in range(self._retry):
                try:
                    cf = ac.run(self.model_name, U=u, runner=runner, **continuation_kwargs)
                except AUTORuntimeError:
                    print(traceback.format_exc())
                    warnings.warn('AUTO continuation failed, possibly retrying.')
                else:
                    break
            else:
                warnings.warn('Problem to complete the forward AUTO continuation, returning nothing.')
                cf = None

        if not self.continuation:
            self.continuation['backward'] = None

        self.continuation['forward'] = cf

        if self.branch_number is None:
            self.branch_number = abs(self.continuation['forward'].data[0].BR)

        if auto_suffix:
            self.auto_save(auto_suffix)

    def make_backward_continuation(self, initial_data, auto_suffix="", **continuation_kwargs):
        runner = ra.runAUTO()
        ac.load(self.model_name, runner=runner)

        if 'MXBF' in continuation_kwargs:
            warnings.warn('Disabling automatic continuation of branch points (MXBF set to 0)')
        continuation_kwargs['MXBF'] = 0

        if 'IBR' not in continuation_kwargs and self.branch_number is not None:
            continuation_kwargs['IBR'] = self.branch_number

        self.initial_data = initial_data

        if isinstance(initial_data, AUTOSolution):
            if 'DS' in continuation_kwargs:
                continuation_kwargs['DS'] = - continuation_kwargs['DS']
                for retry in range(self._retry):
                    try:
                        cb = ac.run(initial_data, runner=runner, **continuation_kwargs)
                    except AUTORuntimeError:
                        print(traceback.format_exc())
                        warnings.warn('AUTO continuation failed, possibly retrying.')
                    else:
                        break
                else:
                    warnings.warn('Problem to complete the backward AUTO continuation, returning nothing.')
                    cb = None
            else:
                for retry in range(self._retry):
                    try:
                        cb = ac.run(initial_data, DS='-', runner=runner, **continuation_kwargs)
                    except AUTORuntimeError:
                        print(traceback.format_exc())
                        warnings.warn('AUTO continuation failed, possibly retrying.')
                    else:
                        break
                else:
                    warnings.warn('Problem to complete the backward AUTO continuation, returning nothing.')
                    cb = None
        else:
            u = {i + 1: initial_data[i] for i in range(self.config_object.ndim)}
            if 'DS' in continuation_kwargs:
                continuation_kwargs['DS'] = - continuation_kwargs['DS']
                for retry in range(self._retry):
                    try:
                        cb = ac.run(self.model_name, U=u, runner=runner, **continuation_kwargs)
                    except AUTORuntimeError:
                        print(traceback.format_exc())
                        warnings.warn('AUTO continuation failed, possibly retrying.')
                    else:
                        break
                else:
                    warnings.warn('Problem to complete the backward AUTO continuation, returning nothing.')
                    cb = None
            else:
                for retry in range(self._retry):
                    try:
                        cb = ac.run(self.model_name, DS='-', U=u, runner=runner, **continuation_kwargs)
                    except AUTORuntimeError:
                        print(traceback.format_exc())
                        warnings.warn('AUTO continuation failed, possibly retrying.')
                    else:
                        break
                else:
                    warnings.warn('Problem to complete the backward AUTO continuation, returning nothing.')
                    cb = None

        if not self.continuation:
            self.continuation['forward'] = None

        self.continuation['backward'] = cb

        if self.branch_number is None:
            self.branch_number = abs(self.continuation['backward'].data[0].BR)

        if auto_suffix:
            self.auto_save(auto_suffix)

    def point_stability(self, idx):
        if isinstance(idx, str):
            if idx[0] == '-':
                if self.continuation['backward'] is not None:
                    s = self.get_solution_by_label(idx)
                    idx = s['PT']
                    ix_map = self._solutions_index_map(direction='backward')
                    if idx is not None:
                        return self.continuation['backward'].data[0].diagnostics[ix_map[idx]]['Eigenvalues']
                    else:
                        warnings.warn('No point stability to show.')
                        return None
                else:
                    warnings.warn('No backward branch to show the stability for.')
                    return None
            else:
                if self.continuation['forward'] is not None:
                    s = self.get_solution_by_label(idx)
                    idx = s['PT']
                    ix_map = self._solutions_index_map(direction='forward')
                    if idx is not None:
                        return self.continuation['forward'].data[0].diagnostics[ix_map[idx]]['Eigenvalues']
                    else:
                        warnings.warn('No point stability to show.')
                        return None
                else:
                    warnings.warn('No forward branch to show the stability for.')
                    return None

        if idx >= 0:
            if self.continuation['forward'] is not None:
                ix_map = self._solutions_index_map(direction='forward')
                if idx in ix_map:
                    return self.continuation['forward'].data[0].diagnostics[ix_map[idx]]['Eigenvalues']
                else:
                    warnings.warn('Point index not found. No point stability to show.')
                    return None
            else:
                warnings.warn('No forward branch to show the stability for.')
                return None
        else:
            if self.continuation['backward'] is not None:
                ix_map = self._solutions_index_map(direction='backward')
                if -idx in ix_map:
                    return self.continuation['backward'].data[0].diagnostics[ix_map[-idx]]['Eigenvalues']
                else:
                    warnings.warn('Point index not found. No point stability to show.')
                    return None
            else:
                warnings.warn('No backward branch to show the stability for.')
                return None

    def _set_from_dict(self, state, load_initial_data=True):
        # store the pathname to pass to the updated class
        state['_path_name'] = self._path_name

        self.__dict__.clear()
        self.__dict__.update(state)
        if isinstance(self.initial_data, dict) and load_initial_data:
            branch_number = abs(self.initial_data['BR'])
            fp_file_list = glob.glob('fp*.pickle')
            fp_branch_numbers = list(map(lambda filename: int(filename.split('_')[1].split('.')[0]), fp_file_list))
            if branch_number in fp_branch_numbers:
                fp = FixedPointContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
                try:
                    fp.load('fp_'+str(branch_number)+'.pickle', load_initial_data=False)

                    for s in fp.full_solutions_list:
                        if s['PT'] == self.initial_data['PT']:
                            self.initial_data = s
                            break
                except FileNotFoundError:
                    warnings.warn('Unable to load initial data. Parent branch was not saved.')
                    self.initial_data = None
            else:
                po_file_list = glob.glob('po*.pickle')
                po_branch_numbers = list(map(lambda s: int(s.split('_')[1].split('.')[0]), po_file_list))
                if branch_number in po_branch_numbers:
                    hp = poc.PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
                    try:
                        hp.load('po_' + str(branch_number) + '.pickle', load_initial_data=False)

                        for s in hp.full_solutions_list:
                            if s['PT'] == self.initial_data['PT']:
                                self.initial_data = s
                                break
                    except FileNotFoundError:
                        warnings.warn('Unable to load initial data. Parent branch was not saved.')
                        self.initial_data = None
                else:
                    warnings.warn('Unable to find initial data.')
                    self.initial_data = None

        if self.auto_filename_suffix:
            self.auto_load(self.auto_filename_suffix)
        else:
            warnings.warn('No AUTO filename suffix specified. Unable to load data.')

    @property
    def isfixedpoint(self):
        return True

    @property
    def isperiodicorbit(self):
        return False
