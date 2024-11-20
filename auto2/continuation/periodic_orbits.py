import os
import sys
import warnings
import logging
import traceback
import glob
from copy import deepcopy

logger = logging.getLogger('logger')

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
    logger.warning('Unable to find auto directory environment variable.')

import auto.AUTOCommands as ac
import auto.runAUTO as ra
from auto.parseS import AUTOSolution
from auto.AUTOExceptions import AUTORuntimeError
from auto2.continuation.base import Continuation

import auto2.continuation.fixed_points as fpc


class PeriodicOrbitContinuation(Continuation):

    def __init__(self, model_name, config_object):

        Continuation.__init__(self, model_name, config_object)

        # plots default behaviours
        self._default_marker = ''
        self._default_markersize = 6.
        self._default_linestyle = '-'
        self._default_linewidth = 1.2

    def make_continuation(self, initial_data, auto_suffix="", only_forward=True, max_bp=None, **continuation_kwargs):

        self.make_forward_continuation(initial_data, "", max_bp=max_bp, **continuation_kwargs)
        if not only_forward:
            self.make_backward_continuation(initial_data, "", max_bp=max_bp, **continuation_kwargs)

        if auto_suffix:
            self.auto_save(auto_suffix)

    def make_forward_continuation(self, initial_data, auto_suffix="", max_bp=None, **continuation_kwargs):
        runner = ra.runAUTO()
        ac.load(self.model_name, runner=runner)

        if 'MXBF' in continuation_kwargs:
            warnings.warn('Disabling automatic continuation of branch points (MXBF set to 0)')
        continuation_kwargs['MXBF'] = 0

        if 'IBR' not in continuation_kwargs and self.branch_number is not None:
            continuation_kwargs['IBR'] = self.branch_number

        if max_bp is not None:
            warnings.warn('Disabling branching points detection after ' + str(max_bp) + ' branching points.')
            if 'SP' in continuation_kwargs:
                continuation_kwargs['SP'].append('BP' + str(max_bp))
            else:
                continuation_kwargs['SP'] = ['BP'+str(max_bp)]

        self.initial_data = initial_data

        if isinstance(initial_data, AUTOSolution):
            for retry in range(self._retry):
                try:
                    cf = ac.run(initial_data, runner=runner, **continuation_kwargs)
                    if max_bp is not None and cf.getIndex(-1)['TY name'] == 'BP':
                        recontinuation_kwargs = deepcopy(continuation_kwargs)
                        for i, sp in enumerate(recontinuation_kwargs['SP']):
                            if 'BP' in sp:
                                recontinuation_kwargs['SP'].pop(i)
                        recontinuation_kwargs['SP'].append('BP0')
                        recontinuation_kwargs['IRS'] = 'BP' + str(max_bp)
                        recontinuation_kwargs['ISW'] = 1
                        recontinuation_kwargs['LAB'] = cf.getIndex(-1)['LAB'] + 1
                        cf2 = ac.run(runner=runner, **recontinuation_kwargs)
                        cf.data[0].append(cf2.data[0])
                except AUTORuntimeError:
                    print(traceback.format_exc())
                    warnings.warn('AUTO continuation failed, possibly retrying.')
                else:
                    break
            else:
                warnings.warn('Problem to complete the forward AUTO continuation, returning nothing.')
                cf = None

        else:
            for retry in range(self._retry):
                try:
                    cf = ac.run(self.model_name, dat=initial_data, runner=runner, **continuation_kwargs)
                    if max_bp is not None and cf.getIndex(-1)['TY name'] == 'BP':
                        recontinuation_kwargs = deepcopy(continuation_kwargs)
                        for i, sp in enumerate(recontinuation_kwargs['SP']):
                            if 'BP' in sp:
                                recontinuation_kwargs['SP'].pop(i)
                        recontinuation_kwargs['SP'].append('BP0')
                        recontinuation_kwargs['SP'].append('BP0')
                        recontinuation_kwargs['IRS'] = 'BP' + str(max_bp)
                        recontinuation_kwargs['ISW'] = 1
                        recontinuation_kwargs['LAB'] = cf.getIndex(-1)['LAB'] + 1
                        cf2 = ac.run(runner=runner, **recontinuation_kwargs)
                        cf.data[0].append(cf2.data[0])
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

    def make_backward_continuation(self, initial_data, auto_suffix="", max_bp=None, **continuation_kwargs):
        runner = ra.runAUTO()
        ac.load(self.model_name, runner=runner)

        if 'MXBF' in continuation_kwargs:
            warnings.warn('Disabling automatic continuation of branch points (MXBF set to 0)')
        continuation_kwargs['MXBF'] = 0

        if 'IBR' not in continuation_kwargs and self.branch_number is not None:
            continuation_kwargs['IBR'] = self.branch_number

        if max_bp is not None:
            warnings.warn('Disabling branching points detection after ' + str(max_bp) + ' branching points.')
            if 'SP' in continuation_kwargs:
                continuation_kwargs['SP'].append('BP' + str(max_bp))
            else:
                continuation_kwargs['SP'] = ['BP'+str(max_bp)]

        self.initial_data = initial_data

        if isinstance(initial_data, AUTOSolution):

            for retry in range(self._retry):
                try:
                    if 'DS' in continuation_kwargs:
                        continuation_kwargs['DS'] = - continuation_kwargs['DS']
                        cb = ac.run(initial_data, runner=runner, **continuation_kwargs)
                    else:
                        cb = ac.run(initial_data, DS='-', runner=runner, **continuation_kwargs)

                    if max_bp is not None and cb.getIndex(-1)['TY name'] == 'BP':
                        recontinuation_kwargs = deepcopy(continuation_kwargs)
                        for i, sp in enumerate(recontinuation_kwargs['SP']):
                            if 'BP' in sp:
                                recontinuation_kwargs['SP'].pop(i)
                        recontinuation_kwargs['SP'].append('BP0')
                        recontinuation_kwargs['IRS'] = 'BP' + str(max_bp)
                        recontinuation_kwargs['ISW'] = 1
                        recontinuation_kwargs['LAB'] = cb.getIndex(-1)['LAB'] + 1
                        cb2 = ac.run(runner=runner, **recontinuation_kwargs)
                        cb.data[0].append(cb2.data[0])
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
                    if 'DS' in continuation_kwargs:
                        continuation_kwargs['DS'] = - continuation_kwargs['DS']
                        cb = ac.run(self.model_name, dat=initial_data, runner=runner, **continuation_kwargs)
                    else:
                        cb = ac.run(self.model_name, DS='-', dat=initial_data, runner=runner, **continuation_kwargs)

                    if max_bp is not None and cb.getIndex(-1)['TY name'] == 'BP':
                        recontinuation_kwargs = deepcopy(continuation_kwargs)
                        for i, sp in enumerate(recontinuation_kwargs['SP']):
                            if 'BP' in sp:
                                recontinuation_kwargs['SP'].pop(i)
                        recontinuation_kwargs['SP'].append('BP0')
                        recontinuation_kwargs['IRS'] = 'BP' + str(max_bp)
                        recontinuation_kwargs['ISW'] = 1
                        recontinuation_kwargs['LAB'] = cb.getIndex(-1)['LAB'] + 1
                        cb2 = ac.run(runner=runner, **recontinuation_kwargs)
                        cb.data[0].append(cb2.data[0])
                except AUTORuntimeError:
                    print(traceback.format_exc())
                    warnings.warn('AUTO continuation failed, possibly retrying.')
                else:
                    break
            else:
                warnings.warn('Problem to complete the backward AUTO continuation, returning only a part.')
                cb = None

        if not self.continuation:
            self.continuation['forward'] = None

        self.continuation['backward'] = cb

        if self.branch_number is None:
            self.branch_number = abs(self.continuation['backward'].data[0].BR)

        if auto_suffix:
            self.auto_save(auto_suffix)

    def orbit_stability(self, idx):
        if isinstance(idx, str):
            if idx[0] == '-':
                if self.continuation['backward'] is not None:
                    s = self.get_solution_by_label(idx)
                    idx = s['PT']
                    ix_map = self._solutions_index_map(direction='backward')
                    if idx is not None:
                        return self.continuation['backward'].data[0].diagnostics[ix_map[idx]]['Multipliers']
                    else:
                        warnings.warn('No orbit stability to show.')
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
                        return self.continuation['forward'].data[0].diagnostics[ix_map[idx]]['Multipliers']
                    else:
                        warnings.warn('No orbit stability to show.')
                        return None
                else:
                    warnings.warn('No forward branch to show the stability for.')
                    return None

        if idx >= 0:
            if self.continuation['forward'] is not None:
                ix_map = self._solutions_index_map(direction='forward')
                if idx in ix_map:
                    return self.continuation['forward'].data[0].diagnostics[ix_map[idx]]['Multipliers']
                else:
                    warnings.warn('Point index not found. No orbit stability to show.')
                    return None
            else:
                warnings.warn('No forward branch to show the stability for.')
                return None
        else:
            if self.continuation['backward'] is not None:
                ix_map = self._solutions_index_map(direction='backward')
                if -idx in ix_map:
                    return self.continuation['backward'].data[0].diagnostics[ix_map[-idx]]['Multipliers']
                else:
                    warnings.warn('Point index not found. No orbit stability to show.')
                    return None
            else:
                warnings.warn('No backward branch to show the stability for.')
                return None

    def _set_from_dict(self, state, load_initial_data=True):
        self.__dict__.clear()
        self.__dict__.update(state)
        if isinstance(self.initial_data, dict) and load_initial_data:
            branch_number = abs(self.initial_data['BR'])
            fp_file_list = glob.glob('fp*.pickle')
            fp_branch_numbers = list(map(lambda filename: int(filename.split('_')[1].split('.')[0]), fp_file_list))
            if branch_number in fp_branch_numbers:
                fp = fpc.FixedPointContinuation(self.model_name, self.config_object)
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
                po_branch_numbers = list(map(lambda filename: int(filename.split('_')[1].split('.')[0]), po_file_list))
                if branch_number in po_branch_numbers:
                    hp = PeriodicOrbitContinuation(self.model_name, self.config_object)
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
        return False

    @property
    def isperiodicorbit(self):
        return True
