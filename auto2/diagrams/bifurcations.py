import os
import sys
import warnings
import logging
import pickle
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.colors import TABLEAU_COLORS

from auto2.parsers.config import ConfigParser
logger = logging.getLogger('logger')
logger.info('Using auto-AUTO (AUTO² or auto2) -- An AUTO-07p automatic search algorithm codebase')
logger.info('Read AUTO-07p manual first before using it. Wishing you a happy continuation, have fun !')
logger.info('Logging messages can be found in the file "auto2.log"')

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

from auto2.continuation.fixed_points import FixedPointContinuation
from auto2.continuation.periodic_orbits import PeriodicOrbitContinuation
from auto.parseS import AUTOSolution


# TODO:  - Add PO computation restart at a given level

class BifurcationDiagram(object):

    def __init__(self, model_name=None):

        self.initial_points = None
        self.model_name = model_name
        if model_name is not None:
            self.config_object = ConfigParser('c.'+model_name)
        else:
            self.config_object = None
        self.fp_branches = dict()
        self.po_branches = dict()

        self.fp_parent = dict()
        self.fp_branches_with_all_bp_computed = list()
        self.fp_branches_with_all_hb_computed = list()
        self.computed_bp_by_fp_branch = dict()
        self.po_parent = dict()
        self.po_branches_with_all_bp_computed = list()
        self.po_branches_with_all_pd_computed = list()
        self.computed_bp_by_po_branch = dict()
        self.computed_pd_by_po_branch = dict()

        self.fp_computed = False
        self.po_computed = False  # maybe not needed

        self.level_reached = 0

        self._comparison_solutions_types = ('HB', 'BP', 'UZ', 'PD', 'EP', 'TR', 'LP')
        self._figure_legend_handles = list()
        self._figure_3d_legend_handles = list()

    def compute_fixed_points_diagram(self, initial_points=None, extra_comparison_parameters=None, comparison_tol=2.e-2, **continuation_kwargs):

        logger.info('Starting the computation of the fixed points bifurcation diagram with model '+str(self.model_name))

        if self.fp_computed:
            logger.warning('Fixed point bifurcation diagram already computed. Aborting.')
            return None

        if initial_points is not None:
            self.initial_points = initial_points

        if 'MNX' not in continuation_kwargs:
            continuation_kwargs['NMX'] = 9000
            logger.info('NMX parameters was not set, so setting it to 9000 points.')

        logger.info('Continuing provided fixed points.')
        br_num = 1
        ncomp = 1
        for point in initial_points:
            parameters = point['parameters']
            initial_data = point['initial_data']

            used_continuation_kwargs = deepcopy(continuation_kwargs)
            if 'PAR' not in used_continuation_kwargs:
                used_continuation_kwargs['PAR'] = dict()
            for par, val in parameters.items():
                if par not in used_continuation_kwargs['PAR']:
                    used_continuation_kwargs['PAR'][par] = val

            used_continuation_kwargs['IBR'] = br_num
            if 'ICP' not in used_continuation_kwargs:
                used_continuation_kwargs['ICP'] = [self.config_object.continuation_parameters[0]]

            fp = FixedPointContinuation(model_name=self.model_name, config_object=self.config_object)
            fp.make_continuation(initial_data, **used_continuation_kwargs)

            self._check_fp_continuation_against_itself(ncomp, fp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

            valid_branch = self._check_fp_continuation_against_other_fp_branches(ncomp, fp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

            if valid_branch:
                self.fp_branches[fp.branch_number] = {'parameters': parameters, 'continuation': fp, 'continuation_kwargs': used_continuation_kwargs}
                self.fp_parent[fp.branch_number] = None
                logger.info('Saving valid branch ' + str(br_num) + ' emanating from detected fixed point ' + str(ncomp) + '.')
                br_num += 1
            ncomp += 1

        logger.info('Provided fixed points continuation has ended.')
        logger.info('Eventually continuing detected branching points.')

        while True:
            nrecomp = 0

            new_branches = dict()

            for parent_branch_number, branch in self.fp_branches.items():

                if parent_branch_number not in self.fp_branches_with_all_bp_computed:

                    forward_branching_points = branch['continuation'].get_filtered_solutions_list(labels='BP', forward=True)
                    backward_branching_points = branch['continuation'].get_filtered_solutions_list(labels='BP', forward=False)

                    direction_and_branching_points = ((1, forward_branching_points), (-1, backward_branching_points))

                    for direction, branching_points in direction_and_branching_points:

                        for ibp, bp in enumerate(branching_points):

                            used_continuation_kwargs = deepcopy(continuation_kwargs)
                            used_continuation_kwargs['ISW'] = -1
                            used_continuation_kwargs['PAR'] = {}

                            for bn in self.computed_bp_by_fp_branch:
                                found_solution = False
                                for ibpt in self.computed_bp_by_fp_branch[bn]:
                                    if ibpt >= 0:
                                        bpt = self.get_continuation(bn).get_solution_by_label('BP' + str(ibpt))
                                    else:
                                        bpt = self.get_continuation(bn).get_solution_by_label('-BP' + str(-ibpt))
                                    if self._check_if_solutions_are_close(bp, bpt, extra_comparison_parameters, comparison_tol):
                                        found_solution = True
                                        break
                                if found_solution:
                                    break
                            else:
                                used_continuation_kwargs['IBR'] = br_num
                                fp = FixedPointContinuation(model_name=self.model_name, config_object=self.config_object)
                                fp.make_continuation(bp, **used_continuation_kwargs)

                                self._check_fp_continuation_against_itself(ncomp, fp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)
                                valid_branch = self._check_fp_continuation_against_other_fp_branches(ncomp, fp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                                if valid_branch:
                                    new_branches[abs(fp.branch_number)] = {'parameters': bp.PAR, 'continuation': fp, 'continuation_kwargs': used_continuation_kwargs}
                                    self.fp_parent[abs(fp.branch_number)] = parent_branch_number
                                    logger.info('Saving valid branch ' + str(br_num) + ' emanating from branch ' + str(parent_branch_number) + '.')
                                    br_num += 1
                                if parent_branch_number not in self.computed_bp_by_fp_branch:
                                    self.computed_bp_by_fp_branch[parent_branch_number] = list()
                                self.computed_bp_by_fp_branch[parent_branch_number].append(direction * (ibp+1))
                                ncomp += 1

                    self.fp_branches_with_all_bp_computed.append(parent_branch_number)
                    nrecomp += 1

            self.fp_branches.update(new_branches)

            if nrecomp == 0:
                break

        logger.info('Fixed points bifurcation diagram computation is over.')
        self.fp_computed = True

    def compute_periodic_orbits_diagram(self, end_level=None, extra_comparison_parameters=None, comparison_tol=2.e-2,
                                        remove_dubious_bp=True, max_number_bp=None, backward_bp_continuation=False, **continuation_kwargs):

        logger.info('Starting the computation of the periodic orbits bifurcation diagram with model '+str(self.model_name))

        if not self.fp_computed:
            logger.warning('Fixed points diagram not computed. No initial data to start with.\n'
                           'Aborting...')
            return None

        if 'MNX' not in continuation_kwargs:
            continuation_kwargs['NMX'] = 9000
            logger.info('NMX parameters was not set, so setting it to 9000 points.')

        br_num = max(self.fp_branches.keys()) + 1
        self.level_reached = 0

        logger.info('First, beginning computation of the periodic orbits from Hopf bifurcations.')
        for parent_branch_number, branch in self.fp_branches.items():

            hb_list = branch['continuation'].get_filtered_solutions_list(labels='HB')

            for hb in hb_list:
                used_continuation_kwargs = deepcopy(continuation_kwargs)
                used_continuation_kwargs['IBR'] = br_num

                if 'IPS' not in used_continuation_kwargs:
                    used_continuation_kwargs['IPS'] = 2

                if 'ICP' not in used_continuation_kwargs:
                    if 11 in self.config_object.parameters_dict.keys():
                        used_continuation_kwargs['ICP'] = [self.config_object.continuation_parameters[0], self.config_object.parameters_dict[11]]
                    else:
                        used_continuation_kwargs['ICP'] = [self.config_object.continuation_parameters[0], 11]
                else:
                    if 11 in self.config_object.parameters_dict.keys():
                        if self.config_object.parameters_dict[11] not in used_continuation_kwargs['ICP'] or 11 not in used_continuation_kwargs['ICP']:
                            used_continuation_kwargs['ICP'].append(self.config_object.parameters_dict[11])
                    else:
                        if 11 not in used_continuation_kwargs['ICP']:
                            used_continuation_kwargs['ICP'].append(11)

                hp = PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object)
                hp.make_continuation(hb, **used_continuation_kwargs)

                self._check_po_continuation_against_itself(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                self._check_po_continuation_against_other_fp_branches(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)
                valid_branch = self._check_po_continuation_against_other_po_branches(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                if valid_branch:
                    self.po_branches[abs(hp.branch_number)] = {'parameters': hb.PAR, 'continuation': hp, 'continuation_kwargs': used_continuation_kwargs}
                    self.po_parent[abs(hp.branch_number)] = parent_branch_number
                    logger.info('Saving valid branch ' + str(br_num) + ' emanating from branch ' + str(parent_branch_number) + '.')
                    br_num += 1
            self.fp_branches_with_all_hb_computed.append(parent_branch_number)

        logger.info('Continuation of the periodic orbits from Hopf bifurcations has ended.')
        self.level_reached += 1
        if self.level_reached == end_level:
            logger.info('As demanded, finishing computation at level '+str(self.level_reached)+' ...')
            return

        logger.info('Beginning computation of the periodic orbits from detected branching and period doubling points.')

        while True:
            nrecomp = 0

            new_branches = dict()

            logger.info('Entering level ' + str(self.level_reached + 1) + ' of continuation...')

            for parent_branch_number, branch in self.po_branches.items():
                logger.debug('Continuing branching points of branch: ' + str(parent_branch_number))

                if parent_branch_number not in self.po_branches_with_all_bp_computed:

                    forward_branching_points = branch['continuation'].get_filtered_solutions_list(labels='BP', forward=True)
                    backward_branching_points = branch['continuation'].get_filtered_solutions_list(labels='BP', forward=False)

                    if max_number_bp is not None:
                        forward_branching_points = forward_branching_points[:max_number_bp]
                        backward_branching_points = backward_branching_points[:max_number_bp]

                    direction_and_branching_points = ((1, forward_branching_points), (-1, backward_branching_points))

                    for direction, branching_points in direction_and_branching_points:

                        for ibp, bp in enumerate(branching_points):

                            used_continuation_kwargs = deepcopy(self.po_branches[parent_branch_number]['continuation_kwargs'])
                            used_continuation_kwargs['ISW'] = -1
                            used_continuation_kwargs['NMX'] = continuation_kwargs['NMX']

                            for bn in self.computed_bp_by_po_branch:
                                found_solution = False
                                for ibpt in self.computed_bp_by_po_branch[bn]:
                                    if ibpt >= 0:
                                        bpt = self.get_continuation(bn).get_solution_by_label('BP' + str(ibpt))
                                    else:
                                        bpt = self.get_continuation(bn).get_solution_by_label('-BP' + str(-ibpt))
                                    if self._check_if_solutions_are_close(bp, bpt, extra_comparison_parameters, comparison_tol):
                                        found_solution = True
                                        break
                                if found_solution:
                                    break
                            else:
                                parent_continuation = self.get_continuation(parent_branch_number)
                                try:
                                    bp_stability = np.array(parent_continuation.orbit_stability(direction * (bp['PT'] - 1)))
                                    # max_accept = 1. / np.nanmin(np.abs(np.where(bp_stability == 0, np.nan, bp_stability)))
                                    max_accept = np.finfo(np.float64).max * 1.e-10
                                    looks_dubious = np.max(bp_stability) > max_accept
                                except ValueError:
                                    par_lst = parent_continuation.continuation_parameters
                                    par_val = [bp.PAR[p] for p in parent_continuation.continuation_parameters]
                                    ini_msg = str(par_lst) + " = " + str(par_val) + ' (branch '+str(parent_branch_number) + ')'
                                    logger.error('No stability information found for PO point at ' + ini_msg + '. Something is wrong, not doing the continuation.')
                                    looks_dubious = True

                                if looks_dubious and remove_dubious_bp:
                                    par_lst = parent_continuation.continuation_parameters
                                    par_val = [bp.PAR[p] for p in parent_continuation.continuation_parameters]
                                    ini_msg = str(par_lst) + " = " + str(par_val) + ' (branch '+str(parent_branch_number) + ')'
                                    try:
                                        s = str(np.max(bp_stability))
                                    except ValueError:
                                        s = '[ unknown ]'
                                    logger.info('Not saving results of PO point at ' + ini_msg + ' because it looks dubious. (max Floquet: ' + s + ' ).'
                                                '\nSkipping to next one.')
                                    valid_branch = False
                                else:
                                    used_continuation_kwargs['IBR'] = br_num
                                    hp = PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object)
                                    hp.make_continuation(bp, only_forward=not backward_bp_continuation, IBR=br_num, ISW=-1)
                                    self._check_po_continuation_against_itself(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                                    self._check_po_continuation_against_other_fp_branches(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)
                                    valid_branch = self._check_po_continuation_against_other_po_branches(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                                if valid_branch:
                                    new_branches[abs(hp.branch_number)] = {'parameters': bp.PAR, 'continuation': hp, 'continuation_kwargs': used_continuation_kwargs}
                                    self.po_parent[abs(hp.branch_number)] = parent_branch_number
                                    logger.info('Saving valid branch ' + str(br_num) + ' emanating from branch ' + str(parent_branch_number) + '.')
                                    br_num += 1
                                if parent_branch_number not in self.computed_bp_by_po_branch:
                                    self.computed_bp_by_po_branch[parent_branch_number] = list()
                                self.computed_bp_by_po_branch[parent_branch_number].append(direction * (ibp+1))

                    self.po_branches_with_all_bp_computed.append(parent_branch_number)
                    nrecomp += 1

                logger.debug('Computation of branching points of branch: ' + str(parent_branch_number) + ' has ended.')
                logger.debug('Continuing period doubling points of branch: ' + str(parent_branch_number))

                if parent_branch_number not in self.po_branches_with_all_pd_computed:

                    forward_period_doubling_points = branch['continuation'].get_filtered_solutions_list(labels='PD', forward=True)
                    backward_period_doubling_points = branch['continuation'].get_filtered_solutions_list(labels='PD', forward=False)

                    direction_and_period_doubling_points = ((1, forward_period_doubling_points), (-1, backward_period_doubling_points))

                    for direction, period_doubling_points in direction_and_period_doubling_points:

                        for ipd, pd in enumerate(period_doubling_points):

                            used_continuation_kwargs = deepcopy(self.po_branches[parent_branch_number]['continuation_kwargs'])
                            used_continuation_kwargs['ISW'] = -1
                            used_continuation_kwargs['NMX'] = continuation_kwargs['NMX']

                            for bn in self.computed_pd_by_po_branch:
                                found_solution = False
                                for ipdt in self.computed_pd_by_po_branch[bn]:
                                    if ipdt >= 0:
                                        pdt = self.get_continuation(bn).get_solution_by_label('PD' + str(ipdt))
                                    else:
                                        pdt = self.get_continuation(bn).get_solution_by_label('-PD' + str(-ipdt))
                                    if self._check_if_solutions_are_close(pd, pdt, extra_comparison_parameters, comparison_tol):
                                        found_solution = True
                                        break
                                if found_solution:
                                    break
                            else:
                                used_continuation_kwargs['IBR'] = br_num
                                hp = PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object)
                                hp.make_continuation(pd, IBR=br_num, ISW=-1)
                                self._check_po_continuation_against_itself(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                                self._check_po_continuation_against_other_fp_branches(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)
                                valid_branch = self._check_po_continuation_against_other_po_branches(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                                if valid_branch:
                                    new_branches[abs(hp.branch_number)] = {'parameters': pd.PAR, 'continuation': hp, 'continuation_kwargs': used_continuation_kwargs}
                                    self.po_parent[abs(hp.branch_number)] = parent_branch_number
                                    logger.debug('Saving branch ' + str(br_num) + ' emanating from branch ' + str(parent_branch_number) + '.')
                                    br_num += 1

                                if parent_branch_number not in self.computed_pd_by_po_branch:
                                    self.computed_pd_by_po_branch[parent_branch_number] = list()
                                self.computed_pd_by_po_branch[parent_branch_number].append(direction * (ipd+1))

                    self.po_branches_with_all_pd_computed.append(parent_branch_number)
                    nrecomp += 1

                logger.debug('Computation of period doubling points of branch: ' + str(parent_branch_number) + ' has ended.')
            self.po_branches.update(new_branches)

            if nrecomp == 0:
                break

            self.level_reached += 1
            logger.info('Computation of level ' + str(self.level_reached) + ' of continuation has ended.')
            if self.level_reached == end_level:
                logger.info('As demanded, finishing computation at level ' + str(self.level_reached) + ' ...')
                break

        logger.info('Finished computation at level ' + str(self.level_reached))

    def _get_dict(self):
        state = self.__dict__.copy()
        fp_branches = dict()
        for branch_number in self.fp_branches:
            fp_branches[branch_number] = dict()
            for key in self.fp_branches[branch_number]:
                if key == 'continuation':
                    fp_branches[branch_number][key] = None
                else:
                    fp_branches[branch_number][key] = self.fp_branches[branch_number][key]
        state['fp_branches'] = fp_branches
        po_branches = dict()
        for branch_number in self.po_branches:
            po_branches[branch_number] = dict()
            for key in self.po_branches[branch_number]:
                if key == 'continuation':
                    po_branches[branch_number][key] = None
                else:
                    po_branches[branch_number][key] = self.po_branches[branch_number][key]
        state['po_branches'] = po_branches
        state['config_object'] = None
        return state

    def _set_from_dict(self, state, load_initial_data=True):
        self.__dict__.clear()
        self.__dict__.update(state)
        self.config_object = ConfigParser('c.'+self.model_name)
        for branch_number in self.fp_branches:
            fp = FixedPointContinuation(self.model_name, self.config_object)
            fp.load('fp_' + str(branch_number) + '.pickle', load_initial_data=load_initial_data)
            self.fp_branches[branch_number]['continuation'] = fp
        for branch_number in self.po_branches:
            po = PeriodicOrbitContinuation(self.model_name, self.config_object)
            po.load('po_' + str(branch_number) + '.pickle', load_initial_data=load_initial_data)
            self.po_branches[branch_number]['continuation'] = po

    def _save_fp_branches(self):
        for branch_number in self.fp_branches:
            self.fp_branches[branch_number]['continuation'].save()

    def _save_po_branches(self):
        for branch_number in self.po_branches:
            self.po_branches[branch_number]['continuation'].save()

    def save(self, filename=None, **kwargs):
        self._save_fp_branches()
        self._save_po_branches()
        if filename is None:
            filename = self.model_name + '.pickle'
            warnings.warn('No filename provided. Using a default one: ' + filename)
        state = self._get_dict()
        with open(filename, 'wb') as f:
            pickle.dump(state, f, **kwargs)

    def load(self, filename, load_initial_data=True, **kwargs):
        try:
            with open(filename, 'rb') as f:
                tmp_dict = pickle.load(f, **kwargs)
        except FileNotFoundError:
            warnings.warn('File not found. Unable to load data.')
            return None

        self._set_from_dict(tmp_dict, load_initial_data=load_initial_data)

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
        for i, param in enumerate(comparison_parameters):
            if isinstance(sol1[param], np.ndarray):
                ssol1[i] = max(sol1[param])
            else:
                ssol1[i] = float(sol1[param])
            if isinstance(sol2[param], np.ndarray):
                ssol2[i] = max(sol2[param])
            else:
                ssol2[i] = float(sol2[param])
        diff = ssol1 - ssol2
        return np.all(np.abs(diff) < tol)

    def _check_fp_continuation_against_itself(self, ncomp, continuation, continuation_kwargs, extra_comparison_parameters, tol):

        fp = continuation
        initial_data = fp.initial_data

        cpar = continuation_kwargs['ICP'][0]
        if isinstance(cpar, int) and self.config_object.parameters_dict:
            cpar = self.config_object.parameters_dict[cpar]
        if extra_comparison_parameters is not None:
            cpar_list = [cpar]
            cpar_list.extend(extra_comparison_parameters)
        else:
            cpar_list = [cpar]

        if fp.continuation['forward'] is not None:

            non_repeating, repeating_solutions = fp.check_for_repetitions(cpar_list, tol=tol, return_repeating_solutions=True, forward=True)
            recompute = False
            for i, v in enumerate(non_repeating[:-1]):
                if not v and not non_repeating[i+1]:
                    recompute = True
                    break

            if recompute:

                first_repeating_sol = repeating_solutions[-1]
                nmx = first_repeating_sol['PT'] + 1
                continuation_kwargs['NMX'] = nmx
                logger.info('Not storing full results of initial point ' + str(ncomp) + ' because it repeats itself (forward).'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                fp.make_forward_continuation(initial_data, **continuation_kwargs)

        if fp.continuation['backward'] is not None:

            non_repeating, repeating_solutions = fp.check_for_repetitions(cpar_list, tol=tol, return_repeating_solutions=True, forward=False)
            recompute = False
            for i, v in enumerate(non_repeating[:-1]):
                if not v and not non_repeating[i+1]:
                    recompute = True
                    break

            if recompute:

                first_repeating_sol = repeating_solutions[-1]
                continuation_kwargs['NMX'] = first_repeating_sol['PT'] + 1
                nmx = first_repeating_sol['PT'] + 1
                continuation_kwargs['NMX'] = nmx
                logger.info('Not storing full results of initial point ' + str(ncomp) + ' because it repeats itself (forward).'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                fp.make_backward_continuation(initial_data, **continuation_kwargs)

    def _check_fp_continuation_against_other_fp_branches(self, ncomp, continuation, continuation_kwargs, extra_comparison_parameters, tol):

        fp = continuation
        initial_data = fp.initial_data

        valid_branch = True
        for n, psol in self.fp_branches.items():
            cpar = continuation_kwargs['ICP'][0]
            if isinstance(cpar, int) and self.config_object.parameters_dict:
                cpar = self.config_object.parameters_dict[cpar]
            if extra_comparison_parameters is not None:
                cpar_list = [cpar]
                cpar_list.extend(extra_comparison_parameters)
            else:
                cpar_list = [cpar]

            if fp.same_solutions_as(psol['continuation'], cpar_list, tol=tol, solutions_types=self._comparison_solutions_types):
                logger.info('Not saving results of initial point ' + str(ncomp) + ' because it already exists as branch ' + str(n) + '.'
                            '\nSkipping to next one.')
                valid_branch = False
                break
            elif fp.solutions_in(psol['continuation'], cpar_list, tol=tol, solutions_types=self._comparison_solutions_types):
                logger.info('Not saving results of initial point ' + str(ncomp) + ' because it is already in branch ' + str(n) + '.'
                            '\nSkipping to next one.')
                valid_branch = False
                break
            else:
                merge_forward, common_solutions = fp.solutions_part_of(psol['continuation'], cpar_list, tol=tol,
                                                                       return_solutions=True, forward=True, solutions_types=self._comparison_solutions_types)
                if merge_forward:
                    first_sol = common_solutions[0]
                    nmx = first_sol['PT'] + 1
                    if nmx > 2:
                        logger.info('Not storing full results of initial point ' + str(ncomp) + ' because it merges forward with branch ' + str(n) + '.'
                                    '\nSaving only the relevant part. NMX set to ' + str(nmx))
                        continuation_kwargs['NMX'] = nmx
                        fp.make_forward_continuation(initial_data, **continuation_kwargs)
                    else:
                        logger.info('Not saving forward results of initial point ' + str(ncomp) + ' because it is already in branch ' + str(n) + '.')
                        fp.continuation['forward'] = None
                else:
                    cross_forward, sol = fp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol,
                                                                  return_solutions=True, forward=True, solutions_types=self._comparison_solutions_types)
                    if cross_forward:
                        nmx = sol['PT'] + 1
                        if nmx > 2:
                            logger.info('Not storing full results of initial point ' + str(ncomp) + ' because it connects forward to branch ' + str(n) + '.'
                                        '\nSaving only the relevant part. NMX set to ' + str(nmx))
                            continuation_kwargs['NMX'] = nmx
                            fp.make_forward_continuation(initial_data, **continuation_kwargs)
                        else:
                            logger.info('Not saving forward results of initial point ' + str(ncomp) + ' because it is already in branch ' + str(n) + '.')
                            fp.continuation['forward'] = None

                merge_backward, common_solutions = fp.solutions_part_of(psol['continuation'], cpar_list, tol=tol,
                                                                        return_solutions=True, forward=False, solutions_types=self._comparison_solutions_types)
                if merge_backward:
                    first_sol = common_solutions[0]
                    nmx = first_sol['PT'] + 1
                    if nmx > 2:
                        logger.info('Not storing full results of initial point ' + str(ncomp) + ' because it merges backward with branch ' + str(n) + '.'
                                    '\nSaving only the relevant part. NMX set to ' + str(nmx))
                        continuation_kwargs['NMX'] = nmx
                        fp.make_backward_continuation(initial_data, **continuation_kwargs)
                    else:
                        logger.info('Not saving backward results of initial point ' + str(ncomp) + ' because it is already in branch ' + str(n) + '.')
                        fp.continuation['backward'] = None
                else:
                    cross_backward, sol = fp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol,
                                                                   return_solutions=True, forward=False, solutions_types=self._comparison_solutions_types)
                    if cross_backward:
                        nmx = sol['PT'] + 1
                        if nmx > 2:
                            logger.info('Not storing full results of initial point ' + str(ncomp) + ' because it connects backward to branch ' + str(n) + '.'
                                        '\nSaving only the relevant part. NMX set to ' + str(nmx))
                            continuation_kwargs['NMX'] = nmx
                            fp.make_backward_continuation(initial_data, **continuation_kwargs)
                        else:
                            logger.info('Not saving backward results of initial point ' + str(ncomp) + ' because it is already in branch ' + str(n) + '.')
                            fp.continuation['backward'] = None

                if fp.continuation['forward'] is None and fp.continuation['backward'] is None:
                    valid_branch = False

        return valid_branch

    def _check_po_continuation_against_other_fp_branches(self, continuation, continuation_kwargs, extra_comparison_parameters, tol):

        hp = continuation
        initial_data = hp.initial_data

        if isinstance(initial_data, str):
            ini_msg = initial_data
        elif isinstance(initial_data, AUTOSolution):
            par_lst = hp.continuation_parameters
            par_val = [initial_data.PAR[p] for p in continuation.continuation_parameters]
            ini_msg = str(par_lst) + " = " + str(par_val) + ' (branch ' + str(abs(initial_data['BR'])) + ')'
        else:
            ini_msg = '[ unknown ]'

        for n, psol in self.fp_branches.items():
            cpar = continuation_kwargs['ICP'][0]
            if isinstance(cpar, int) and self.config_object.parameters_dict:
                cpar = self.config_object.parameters_dict[cpar]
            if extra_comparison_parameters is not None:
                cpar_list = [cpar]
                cpar_list.extend(extra_comparison_parameters)
            else:
                cpar_list = [cpar]

            crossing, sol = hp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol, return_solutions=True, forward=True,
                                                     solutions_types=self._comparison_solutions_types)
            if crossing and not self._check_if_solutions_are_close(initial_data, sol, extra_comparison_parameters, tol):
                nmx = sol['PT'] + 1
                logger.info('Not storing full results of PO point at ' + ini_msg + ' because it connects to branch ' + str(n) + '.'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                continuation_kwargs['NMX'] = nmx
                hp.make_forward_continuation(initial_data, **continuation_kwargs)

            crossing, sol = hp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol, return_solutions=True, forward=False,
                                                     solutions_types=self._comparison_solutions_types)
            if crossing and not self._check_if_solutions_are_close(initial_data, sol, extra_comparison_parameters, tol):
                nmx = sol['PT'] + 1
                logger.info('Not storing full results of PO point at ' + ini_msg + ' because it connects to branch ' + str(n) + '.'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                continuation_kwargs['NMX'] = nmx
                hp.make_backward_continuation(initial_data, **continuation_kwargs)

    def _check_po_continuation_against_itself(self, continuation, continuation_kwargs, extra_comparison_parameters, tol):

        hp = continuation
        initial_data = hp.initial_data

        if isinstance(initial_data, str):
            ini_msg = initial_data
        elif isinstance(initial_data, AUTOSolution):
            par_lst = hp.continuation_parameters
            par_val = [initial_data.PAR[p] for p in continuation.continuation_parameters]
            ini_msg = str(par_lst) + " = " + str(par_val) + ' (branch ' + str(abs(initial_data['BR'])) + ')'
        else:
            ini_msg = '[ unknown ]'

        cpar = continuation_kwargs['ICP'][0]
        if isinstance(cpar, int) and self.config_object.parameters_dict:
            cpar = self.config_object.parameters_dict[cpar]
        if extra_comparison_parameters is not None:
            cpar_list = [cpar]
            cpar_list.extend(extra_comparison_parameters)
        else:
            cpar_list = [cpar]

        if hp.continuation['forward'] is not None:

            non_repeating, repeating_solutions = hp.check_for_repetitions(cpar_list, tol=tol, return_repeating_solutions=True, forward=True)
            recompute = False
            for i, v in enumerate(non_repeating[:-1]):
                if not v and not non_repeating[i+1]:
                    recompute = True
                    break

            if recompute:
                first_repeating_sol = repeating_solutions[-1]
                nmx = first_repeating_sol['PT'] + 1
                logger.info('Not storing full results of PO point at ' + ini_msg + ' because it repeats itself (forward).'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                continuation_kwargs['NMX'] = nmx
                hp.make_forward_continuation(initial_data, **continuation_kwargs)

        if hp.continuation['backward'] is not None:

            non_repeating, repeating_solutions = hp.check_for_repetitions(cpar_list, tol=tol, return_repeating_solutions=True, forward=False)
            recompute = False
            for i, v in enumerate(non_repeating[:-1]):
                if not v and not non_repeating[i+1]:
                    recompute = True
                    break

            if recompute:
                first_repeating_sol = repeating_solutions[-1]
                nmx = first_repeating_sol['PT'] + 1
                logger.info('Not storing full results of PO point at ' + ini_msg + ' because it repeats itself (backward).'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                continuation_kwargs['NMX'] = nmx
                hp.make_backward_continuation(initial_data, **continuation_kwargs)

    def _check_po_continuation_against_other_po_branches(self, continuation, continuation_kwargs, extra_comparison_parameters, tol):

        hp = continuation
        initial_data = hp.initial_data

        if isinstance(initial_data, str):
            ini_msg = initial_data
        elif isinstance(initial_data, AUTOSolution):
            par_lst = hp.continuation_parameters
            par_val = [initial_data.PAR[p] for p in continuation.continuation_parameters]
            ini_msg = str(par_lst) + " = " + str(par_val) + ' (branch ' + str(abs(initial_data['BR'])) + ')'
        else:
            ini_msg = '[ unknown ]'

        valid_branch = True
        for n, psol in self.po_branches.items():
            cpar = continuation_kwargs['ICP'][0]
            if isinstance(cpar, int) and self.config_object.parameters_dict:
                cpar = self.config_object.parameters_dict[cpar]
            if extra_comparison_parameters is not None:
                cpar_list = [cpar]
                cpar_list.extend(extra_comparison_parameters)
            else:
                cpar_list = [cpar]

            if hp.same_solutions_as(psol['continuation'], cpar_list, tol=tol, solutions_types=self._comparison_solutions_types):
                logger.info('Not saving results of PO point at ' + ini_msg + ' because it already exists (branch ' + str(n) + ').'
                            '\nSkipping to next one.')
                valid_branch = False
                break
            elif hp.solutions_in(psol['continuation'], cpar_list, tol=tol, solutions_types=self._comparison_solutions_types):
                logger.info('Not saving results of PO point at ' + ini_msg + ' because it is already in branch ' + str(n) + '.'
                            '\nSkipping to next one.')
                valid_branch = False
                break
            else:
                merge_forward, common_solutions = hp.solutions_part_of(psol['continuation'], cpar_list, tol=tol,
                                                                       return_solutions=True, forward=True, solutions_types=self._comparison_solutions_types)
                remake_continuation = True
                if merge_forward:
                    first_sol = common_solutions[0]
                    if isinstance(initial_data, AUTOSolution):
                        if psol['continuation'].branch_number == abs(initial_data['BR']) and first_sol['PT'] < 10:
                            remake_continuation = not self._check_if_solutions_are_close(initial_data, first_sol, extra_comparison_parameters, tol)
                    if remake_continuation:
                        nmx = first_sol['PT'] + 1
                        logger.info('Not storing full results of PO point at ' + ini_msg + ' because it merges forward with branch ' + str(n) + '.'
                                    '\nSaving only the relevant part. NMX set to ' + str(nmx))
                        continuation_kwargs['NMX'] = nmx
                        hp.make_forward_continuation(initial_data, **continuation_kwargs)

                if not merge_forward and not remake_continuation:
                    cross_forward, sol = hp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol,
                                                                  return_solutions=True, forward=True, solutions_types=self._comparison_solutions_types)

                    if cross_forward:
                        nmx = sol['PT'] + 1
                        logger.info('Not storing full results of PO point at ' + ini_msg + ' because it connects forward to branch ' + str(n) + '.'
                                    '\nSaving only the relevant part. NMX set to ' + str(nmx))
                        continuation_kwargs['NMX'] = nmx
                        hp.make_forward_continuation(initial_data, **continuation_kwargs)

                merge_backward, common_solutions = hp.solutions_part_of(psol['continuation'], cpar_list, tol=tol,
                                                                        return_solutions=True, forward=False, solutions_types=self._comparison_solutions_types)
                remake_continuation = True
                if merge_backward:
                    first_sol = common_solutions[0]
                    if isinstance(initial_data, AUTOSolution):
                        if psol['continuation'].branch_number == abs(initial_data['BR']) and first_sol['PT'] < 10:
                            remake_continuation = not self._check_if_solutions_are_close(initial_data, first_sol, extra_comparison_parameters, tol)
                    if remake_continuation:
                        nmx = first_sol['PT'] + 1
                        logger.info('Not storing full results of PO point at ' + ini_msg + ' because it merges backward with branch ' + str(n) + '.'
                                    '\nSaving only the relevant part. NMX set to ' + str(nmx))
                        continuation_kwargs['NMX'] = nmx
                        hp.make_backward_continuation(initial_data, **continuation_kwargs)

                if not merge_backward or not remake_continuation:
                    cross_backward, sol = hp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol,
                                                                   return_solutions=True, forward=False, solutions_types=self._comparison_solutions_types)

                    if cross_backward:
                        nmx = sol['PT'] + 1
                        logger.info('Not storing full results of  PO point at ' + ini_msg + ' because it connects backward to branch ' + str(n) + '.'
                                    '\nSaving only the relevant part. NMX set to ' + str(nmx))
                        continuation_kwargs['NMX'] = nmx
                        hp.make_backward_continuation(initial_data, **continuation_kwargs)

        return valid_branch

    def plot_fixed_points_diagram(self, variables=(0, 1), ax=None, figsize=(10, 8), cmap=None, return_used_colors=False,
                                  legend=True, **kwargs):

        if 'plot_kwargs' not in kwargs:
            kwargs['plot_kwargs'] = dict()

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.gca()

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        used_colors = dict()
        colors_list = list(TABLEAU_COLORS.keys())
        new_handles = list()
        for i, b in enumerate(self.fp_branches):
            if cmap is None:
                kwargs['plot_kwargs']['color'] = colors_list[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_fp_branches)
            self.fp_branches[b]['continuation'].plot_branche_parts(variables, ax=ax, **kwargs)
            used_colors[b] = kwargs['plot_kwargs']['color']

        for i, b in enumerate(self.fp_branches):
            if cmap is None:
                kwargs['plot_kwargs']['color'] = colors_list[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_fp_branches)
            kwargs['plot_kwargs']['linestyle'] = '-'
            kwargs['plot_kwargs']['label'] = 'BR '+str(b)
            new_handles.append(Line2D([], [], **(kwargs['plot_kwargs'])))

        if len(self._figure_legend_handles) in (0, self.number_of_po_branches):
            self._figure_legend_handles.extend(new_handles)

        if legend:
            ax.legend(handles=self._figure_legend_handles)

        if return_used_colors:
            return ax, used_colors
        else:
            return ax

    def plot_fixed_points_diagram_3D(self, variables=(0, 1, 2), ax=None, figsize=(10, 8), cmap=None, return_used_colors=False,
                                     legend=True, **kwargs):

        if 'plot_kwargs' not in kwargs:
            kwargs['plot_kwargs'] = dict()

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(projection='3d')

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        used_colors = dict()
        colors_list = list(TABLEAU_COLORS.keys())
        new_handles = list()
        for i, b in enumerate(self.fp_branches):
            if cmap is None:
                kwargs['plot_kwargs']['color'] = colors_list[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_fp_branches)
            self.fp_branches[b]['continuation'].plot_branche_parts_3D(variables, ax=ax, **kwargs)
            used_colors[b] = kwargs['plot_kwargs']['color']

        for i, b in enumerate(self.fp_branches):
            if cmap is None:
                kwargs['plot_kwargs']['color'] = colors_list[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_fp_branches)
            kwargs['plot_kwargs']['linestyle'] = '-'
            kwargs['plot_kwargs']['label'] = 'BR ' + str(b)
            new_handles.append(Line2D([], [], **(kwargs['plot_kwargs'])))

        if len(self._figure_3d_legend_handles) in (0, self.number_of_po_branches):
            self._figure_3d_legend_handles.extend(new_handles)

        if legend:
            ax.legend(handles=self._figure_3d_legend_handles)

        if return_used_colors:
            return ax, used_colors
        else:
            return ax

    def plot_periodic_orbits_diagram(self, variables=(0, 1), ax=None, figsize=(10, 8), cmap=None, return_used_colors=False,
                                     legend=True, **kwargs):

        if 'plot_kwargs' not in kwargs:
            kwargs['plot_kwargs'] = dict()

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.gca()

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        used_colors = dict()
        colors_list = list(TABLEAU_COLORS.keys())
        new_handles = list()
        for i, b in enumerate(self.po_branches):
            if cmap is None:
                kwargs['plot_kwargs']['color'] = colors_list[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_po_branches)
            self.po_branches[b]['continuation'].plot_branche_parts(variables, ax=ax, **kwargs)
            used_colors[b] = kwargs['plot_kwargs']['color']

        for i, b in enumerate(self.po_branches):
            if cmap is None:
                kwargs['plot_kwargs']['color'] = colors_list[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_po_branches)
            kwargs['plot_kwargs']['linestyle'] = '-'
            kwargs['plot_kwargs']['label'] = 'BR '+str(b)
            new_handles.append(Line2D([], [], **(kwargs['plot_kwargs'])))

        if len(self._figure_legend_handles) in (0, self.number_of_fp_branches):
            self._figure_legend_handles.extend(new_handles)

        if legend:
            ax.legend(handles=self._figure_legend_handles)

        if return_used_colors:
            return ax, used_colors
        else:
            return ax

    def plot_periodic_orbits_diagram_3D(self, variables=(0, 1, 2), ax=None, figsize=(10, 8), cmap=None, return_used_colors=False,
                                        legend=True, **kwargs):

        if 'plot_kwargs' not in kwargs:
            kwargs['plot_kwargs'] = dict()

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(projection='3d')

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        used_colors = dict()
        colors_list = list(TABLEAU_COLORS.keys())
        new_handles = list()
        for i, b in enumerate(self.po_branches):
            if cmap is None:
                kwargs['plot_kwargs']['color'] = colors_list[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_po_branches)
            self.po_branches[b]['continuation'].plot_branche_parts_3D(variables, ax=ax, **kwargs)
            used_colors[b] = kwargs['plot_kwargs']['color']

        for i, b in enumerate(self.po_branches):
            if cmap is None:
                kwargs['plot_kwargs']['color'] = colors_list[i]
            else:
                kwargs['plot_kwargs']['color'] = cmap(i / self.number_of_po_branches)
            kwargs['plot_kwargs']['linestyle'] = '-'
            kwargs['plot_kwargs']['label'] = 'BR ' + str(b)
            new_handles.append(Line2D([], [], **(kwargs['plot_kwargs'])))

        if len(self._figure_3d_legend_handles) in (0, self.number_of_fp_branches):
            self._figure_3d_legend_handles.extend(new_handles)

        if legend:
            ax.legend(handles=self._figure_3d_legend_handles)

        if return_used_colors:
            return ax, used_colors
        else:
            return ax

    def plot_diagram_and_solutions(self, solutions_parameter_value, diagram_variables=(1,), solutions_variables=(0, 1),
                                   axes=None, figsize=(10, 16), solutions_tol=0.01, fixed_points_diagram_kwargs=None, periodic_orbits_diagram_kwargs=None,
                                   solutions_kwargs=None):

        if axes is None:
            fig, axes = plt.subplots(2, 1, figsize=figsize)

        if self.fp_branches:
            n = next(iter(self.fp_branches))
            parameter = self.fp_branches[n]['continuation_kwargs']['ICP']
        elif self.po_branches:
            n = next(iter(self.po_branches))
            parameter = self.po_branches[n]['continuation_kwargs']['ICP']
        else:
            return None

        if self.fp_branches:
            if fixed_points_diagram_kwargs is None:
                fixed_points_diagram_kwargs = dict()

            _, fp_used_colors = self.plot_fixed_points_diagram(variables=(parameter, diagram_variables[0]), ax=axes[0],
                                                               return_used_colors=True, **fixed_points_diagram_kwargs)

        if self.po_branches:
            if periodic_orbits_diagram_kwargs is None:
                periodic_orbits_diagram_kwargs = dict()

            _, po_used_colors = self.plot_periodic_orbits_diagram(variables=(parameter, diagram_variables[0]), ax=axes[0],
                                                                  return_used_colors=True, **periodic_orbits_diagram_kwargs)

        axes[0].axvline(x=solutions_parameter_value, linestyle='--', color='k', linewidth=1.2)
        axes[0].set_title('Bifurcation diagram')

        # part on solutions

        if solutions_kwargs is None:
            solutions_kwargs = dict()

        if 'plot_kwargs' not in solutions_kwargs:
            solutions_kwargs['plot_kwargs'] = dict()

        for branch_number in self.fp_branches:

            fp = self.fp_branches[branch_number]['continuation']
            solutions_kwargs['plot_kwargs']['color'] = fp_used_colors[branch_number]
            fp.plot_solutions(solutions_variables, ax=axes[1], parameter=parameter, value=solutions_parameter_value, tol=solutions_tol, **solutions_kwargs)

        for branch_number in self.po_branches:

            hp = self.po_branches[branch_number]['continuation']
            solutions_kwargs['plot_kwargs']['color'] = po_used_colors[branch_number]
            hp.plot_solutions(solutions_variables, ax=axes[1], parameter=parameter, value=solutions_parameter_value, tol=solutions_tol, **solutions_kwargs)

        axes[1].set_title('Solution in phase space')

        return axes

    def plot_diagram_in_3D_and_solutions_in_3D(self, solutions_parameter_value, diagram_variables=(1, 2), solutions_variables=(0, 1, 2),
                                               axes=None, figsize=(10, 16), solutions_tol=0.01, fixed_points_diagram_kwargs=None,
                                               periodic_orbits_diagram_kwargs=None, solutions_kwargs=None):

        if axes is None:
            fig, axes = plt.subplots(2, 1, figsize=figsize, subplot_kw={'projection': '3d'})

        if self.fp_branches:
            n = next(iter(self.fp_branches))
            parameter = self.fp_branches[n]['continuation_kwargs']['ICP']
        elif self.po_branches:
            n = next(iter(self.po_branches))
            parameter = self.po_branches[n]['continuation_kwargs']['ICP']
        else:
            return None

        if self.fp_branches:
            if fixed_points_diagram_kwargs is None:
                fixed_points_diagram_kwargs = dict()

            _, fp_used_colors = self.plot_fixed_points_diagram_3D(variables=(parameter, *diagram_variables), ax=axes[0],
                                                                  return_used_colors=True, **fixed_points_diagram_kwargs)

        if self.po_branches:
            if periodic_orbits_diagram_kwargs is None:
                periodic_orbits_diagram_kwargs = dict()

            _, po_used_colors = self.plot_periodic_orbits_diagram_3D(variables=(parameter, *diagram_variables), ax=axes[0],
                                                                     return_used_colors=True, **periodic_orbits_diagram_kwargs)

        # axes[0].axvline(x=solutions_parameter_value, linestyle='--', color='k', linewidth=1.2)
        axes[0].set_title('Bifurcation diagram')

        # part on solutions

        if solutions_kwargs is None:
            solutions_kwargs = dict()

        if 'plot_kwargs' not in solutions_kwargs:
            solutions_kwargs['plot_kwargs'] = dict()

        for branch_number in self.fp_branches:

            fp = self.fp_branches[branch_number]['continuation']
            solutions_kwargs['plot_kwargs']['color'] = fp_used_colors[branch_number]
            fp.plot_solutions_3D(solutions_variables, ax=axes[1], parameter=parameter, value=solutions_parameter_value, tol=solutions_tol, **solutions_kwargs)

        for branch_number in self.po_branches:

            hp = self.po_branches[branch_number]['continuation']
            solutions_kwargs['plot_kwargs']['color'] = po_used_colors[branch_number]
            hp.plot_solutions_3D(solutions_variables, ax=axes[1], parameter=parameter, value=solutions_parameter_value, tol=solutions_tol, **solutions_kwargs)

        axes[1].set_title('Solution in phase space')

        return axes

    def plot_diagram_and_solutions_in_3D(self, solutions_parameter_value, diagram_variables=(1,), solutions_variables=(0, 1, 2),
                                         axes=None, figsize=(10, 16), solutions_tol=0.01, fixed_points_diagram_kwargs=None,
                                         periodic_orbits_diagram_kwargs=None, solutions_kwargs=None):

        if axes is None:
            axes = list()
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(2, 1, 1)
            axes.append(ax)
            ax = fig.add_subplot(2, 1, 2, projection='3d')
            axes.append(ax)

        if self.fp_branches:
            n = next(iter(self.fp_branches))
            parameter = self.fp_branches[n]['continuation_kwargs']['ICP']
        elif self.po_branches:
            n = next(iter(self.po_branches))
            parameter = self.po_branches[n]['continuation_kwargs']['ICP']
        else:
            return None

        if self.fp_branches:
            if fixed_points_diagram_kwargs is None:
                fixed_points_diagram_kwargs = dict()

            _, fp_used_colors = self.plot_fixed_points_diagram(variables=(parameter, diagram_variables[0]), ax=axes[0],
                                                               return_used_colors=True, **fixed_points_diagram_kwargs)

        if self.po_branches:
            if periodic_orbits_diagram_kwargs is None:
                periodic_orbits_diagram_kwargs = dict()

            _, po_used_colors = self.plot_periodic_orbits_diagram(variables=(parameter, diagram_variables[0]), ax=axes[0],
                                                                  return_used_colors=True, **periodic_orbits_diagram_kwargs)

        axes[0].axvline(x=solutions_parameter_value, linestyle='--', color='k', linewidth=1.2)
        axes[0].set_title('Bifurcation diagram')

        # part on solutions

        if solutions_kwargs is None:
            solutions_kwargs = dict()

        if 'plot_kwargs' not in solutions_kwargs:
            solutions_kwargs['plot_kwargs'] = dict()

        for branch_number in self.fp_branches:

            fp = self.fp_branches[branch_number]['continuation']
            solutions_kwargs['plot_kwargs']['color'] = fp_used_colors[branch_number]
            fp.plot_solutions_3D(solutions_variables, ax=axes[1], parameter=parameter, value=solutions_parameter_value, tol=solutions_tol, **solutions_kwargs)

        for branch_number in self.po_branches:

            hp = self.po_branches[branch_number]['continuation']
            solutions_kwargs['plot_kwargs']['color'] = po_used_colors[branch_number]
            hp.plot_solutions_3D(solutions_variables, ax=axes[1], parameter=parameter, value=solutions_parameter_value, tol=solutions_tol, **solutions_kwargs)

        axes[1].set_title('Solution in phase space')

        return axes

    @property
    def number_of_fp_branches(self):
        return len(self.fp_branches.keys())

    @property
    def number_of_po_branches(self):
        return len(self.po_branches.keys())

    def get_continuation(self, idx):
        if idx in self.fp_branches:
            return self.fp_branches[idx]['continuation']
        elif idx in self.po_branches:
            return self.po_branches[idx]['continuation']
        else:
            return None

    @property
    def periodic_orbits_variables_list(self):
        for branch_number in self.po_branches:
            branch = self.get_continuation(branch_number)
            sl = branch.available_variables
            if sl is not None:
                return deepcopy(sl)
        return list()

    @property
    def fixed_points_variables_list(self):
        for branch_number in self.fp_branches:
            branch = self.get_continuation(branch_number)
            sl = branch.available_variables
            if sl is not None:
                return deepcopy(sl)
        return list()

    @property
    def fixed_points_solutions_list(self):
        sl = list()
        for branch_number in self.fp_branches:
            branch = self.get_continuation(branch_number)
            sl.extend(branch.full_solutions_list)
        return sl

    @property
    def periodic_orbits_solutions_list(self):
        sl = list()
        for branch_number in self.po_branches:
            branch = self.get_continuation(branch_number)
            sl.extend(branch.full_solutions_list)
        return sl

    @property
    def full_solutions_list(self):
        sl = self.fixed_points_solutions_list
        sl.extend(self.periodic_orbits_solutions_list)
        return sl
