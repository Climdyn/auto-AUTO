
"""

Bifurcation diagram class definition
====================================

This module defines the bifurcation diagram object used in auto-AUTO.

"""

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
        sys.path.append(os.path.join(auto_directory, 'python'))
except KeyError:
    logger.warning('Unable to find auto directory environment variable.')

from auto2.continuations.fixed_points import FixedPointContinuation
from auto2.continuations.periodic_orbits import PeriodicOrbitContinuation
from auto.parseS import AUTOSolution


class BifurcationDiagram(object):
    """Base class for any continuation in auto-AUTO.

    Parameters
    ----------
    model_name: str
        The name of the model to load. Used to load the .f90 file of the provided model.
    path_name: str, optional
        The directory path where files are read/saved.
        If `None`, defaults to current working directory.

    Attributes
    ----------
    initial_points: list(dict)
        List of initial points used to start the bifurcation diagram. It has the following structure:
        Each point is an entry in the list and is specified with a dictionary with the following two keys:

        * `'parameters'`: a dictionary of with the used continuation parameter and its value at which the `initial_data` (see below)
          is valid. I.e. its structure is `{'parameter_name': parameter_value}`.
        * `'initial_data'`: a 1-D Numpy :class:`~numpy.ndarray` with the coordinates in phase space of the initial point, or an
          AUTOSolution object.

    model_name: str
        The name of the loaded model. Used to load the .f90 file of the provided model.
    config_object: ~auto2.parsers.config.ConfigParser
        The ConfigParser object used during the continuations.
    fp_branches: dict
        Dictionary holding the continuations data of the fixed points of the bifurcation diagram, each labelled by their
        branch number (i.e. dictionary keys are the branch number).
        Each entry of the dictionary is itself a dictionary characterizing te branch and with the following entries:

        * `'parameters'`: A list of the parameters values used to start the continuation(s).
        * `'continuation'`: A :class:`~auto2.continuations.base.Continuation` object including the continuation(s) data.
        * `'continuation_kwargs'`: A dictionary containing the AUTO parameters used to start the continuation(s).
    po_branches: dict
        Dictionary holding the continuations of the periodic orbits of the bifurcation diagram, each labelled by their
        branch number (i.e. dictionary keys are the branch number).
        Each entry of the dictionary is itself a dictionary characterizing te branch and with the following entries:

        * `'parameters'`: A list of the parameters values used to start the continuation(s).
        * `'continuation'`: A :class:`~auto2.continuations.base.Continuation` object including the continuation(s) data.
        * `'continuation_kwargs'`: A dictionary containing the AUTO parameters used to start the continuation(s).
    fp_parent: dict(int or None)
        Dictionary holding the branch number of the parents of a given fixed point branch (whose branch number is the dictionary key).
        `None` if there is no parent branch.
    fp_branches_with_all_bp_computed: list(int)
        List of fixed point branch number for which all the branching points of the branch have been continued.
    fp_branches_with_all_hb_computed: list(int)
        List of fixed point branch number for which all the Hopf bifurcations of the branch have been continued.
    computed_bp_by_fp_branch: dict(list(int))
        List of solution label number of branching points which have been continued for each fixed point branch
        (indexed by their branch number forming the dictionary keys).
        E.g. for a given branch, `[1, 2]` indicates that the branches emanating from
        the first and second branching points on the `forward` continuation have been computed.
        Negative label numbers indicate branching points on the `backward` continuation.
    computed_hb_by_fp_branch: dict(list(int))
        List of solution label number of Hopf bifurcations which have been continued for each fixed point branch
        (indexed by their branch number forming the dictionary keys).
        E.g. for a given branch, `[1, 2]` indicates that the branches emanating from
        the first and second Hopf points on the `forward` continuation have been computed.
        Negative label numbers indicate branching points on the `backward` continuation.
    po_parent: dict(int or None)
        Dictionary holding the branch number of the parents of a given periodic orbit branch (whose branch number is the dictionary key).
        `None` if there is no parent branch.
    po_branches_with_all_bp_computed: list(int)
        List of periodic orbit branch number for which all the branching points of the branch have been continued.
    po_branches_with_all_pd_computed: list(int)
        List of periodic orbit branch number for which all the period doubling bifurcations of the branch have been continued.
    computed_bp_by_po_branch: dict(list(int))
        List of solution label number of branching points which have been continued for each periodic orbit
        (indexed by their branch number forming the dictionary keys).
        E.g. for a given branch, `[1, 2]` indicates that the branches emanating from
        the first and second branching points on the `forward` continuation have been computed.
        Negative label numbers indicate branching points on the `backward` continuation.
    valid_bp_by_po_branch:
        List of solution label number periodic orbit branching points which have resulted in a valid continuation for each periodic orbit
        (indexed by their branch number forming the dictionary keys).
        E.g. for a given branch, `[1, 2]` indicates that the branches emanating from
        the first and second periodic orbit branching points on the `forward` continuation were valid.
        Negative label numbers indicate branching points on the `backward` continuation.
        **A valid continuation is a continuation for which no subset is identified with another branch or another part of itself (loops).**
    computed_pd_by_po_branch: dict(list(int))
        List of solution label number of period doubling bifurcations which have been continued for each periodic orbit
        (indexed by their branch number forming the dictionary keys).
        E.g. for a given branch, `[1, 2]` indicates that the branches emanating from
        the first and second period doubling points on the `forward` continuation have been computed.
        Negative label numbers indicate branching points on the `backward` continuation.
    valid_pd_by_po_branch: dict(list(int))
        List of solution label number of period doubling bifurcations which have resulted in a valid continuation for each periodic orbit
        (indexed by their branch number forming the dictionary keys).
        E.g. for a given branch, `[1, 2]` indicates that the branches emanating from
        the first and second period doubling points on the `forward` continuation were valid.
        Negative label numbers indicate branching points on the `backward` continuation.
        **A valid continuation is a continuation for which no subset is identified with another branch or another part of itself (loops).**
    fp_computed: bool
        Whether all the fixed points continuations possibilities of the bifurcation diagram have been computed.
    po_computed: bool
        Whether all the periodic orbits continuations possibilities of the bifurcation diagram have been computed.

    """

    def __init__(self, model_name=None, path_name=None):

        self.initial_points = None
        self.model_name = model_name

        self._path_name = None
        if path_name is not None:
            self.set_path_name(path_name)

        if model_name is not None:
            filepath = 'c.'+model_name if self._path_name is None else os.path.join(self._path_name, 'c.'+model_name)
            self.config_object = ConfigParser(filepath)
        else:
            self.config_object = None
            
        self.fp_branches = dict()
        self.po_branches = dict()

        self.fp_parent = dict()
        self.fp_branches_with_all_bp_computed = list()
        self.fp_branches_with_all_hb_computed = list()
        self.computed_bp_by_fp_branch = dict()
        self.computed_hb_by_fp_branch = dict()
        self.po_parent = dict()
        self.po_branches_with_all_bp_computed = list()
        self.po_branches_with_all_pd_computed = list()
        self.computed_bp_by_po_branch = dict()
        self.valid_bp_by_po_branch = dict()
        self.computed_pd_by_po_branch = dict()
        self.valid_pd_by_po_branch = dict()

        self.fp_computed = False
        self.po_computed = False

        self.level_reached = 0

        self._comparison_solutions_types = ('HB', 'BP', 'UZ', 'PD', 'EP', 'TR', 'LP')
        self._figure_legend_handles = list()
        self._figure_3d_legend_handles = list()

    @property
    def path_name(self):
        """str: The path where the |AUTO| and AUTO² files must be or are stored."""
        return self._path_name

    def set_path_name(self, path_name):
        """Set the path where the |AUTO| and AUTO² files must be or are stored.

        Parameters
        ----------
        path_name: str
            The path.
        """
        if os.path.exists(path_name):
            self._path_name = path_name
        else:
            warnings.warn("Path name given does not exist. Using the current working directory.")
            self._path_name = None

    def compute_fixed_points_diagram(self, initial_points=None, extra_comparison_parameters=None, comparison_tol=2.e-2, **continuation_kwargs):
        """Method which starts the computation of a fixed points bifurcation diagrams using a set of provided initial points.

        Parameters
        ----------
        initial_points: list(dict) or None, optional
            List of initial fixed points to consider for computing the diagram.
            Each entry in the list is a dictionary with the following two keys:

            * `'parameters'`: a dictionary of with the used continuation parameter and its value at which the `initial_data` (see below)
              is valid. I.e. its structure is `{'parameter_name': parameter_value}`. Can be overloaded by the values found in
              the `continuation_kwargs`.
            * `'initial_data'`: a 1-D Numpy :class:`~numpy.ndarray` with the coordinates in phase space of the initial point, or an
              AUTOSolution object.

            If `None`, recompute the fixed points already present in the :attr:`initial_points` class attribute.
            Default to `None`.
        extra_comparison_parameters: list(str or int) or None
            List of extra parameters labels or numbers to use in the various comparison to determine the validity of the continuations.
            If `None`, only the defaults parameters will be considered, i.e. the continuation parameter and the `L2` norm.
            Otherwise, the specified parameters are added to the two default ones.
            Default to `None`.
        comparison_tol: float or list(float) or ~numpy.ndarray(float), optional
            The numerical tolerance of the parameters comparison done to check the continuations validity.
            If a single float is provided, assume the same tolerance for each parameter.
            If a list or a 1-D array is passed, it must have the dimension `2+len(extra_comparison_parameters)`, each value in the
            array corresponding to a parameter.
            Default to `0.02`.
        continuation_kwargs: dict
            Parameters to pass to AUTO for each continuation.
            See the :meth:`~auto2.continuations.fixed_points.FixedPointContinuation.make_continuation` method of the
            :class:`~~auto2.continuations.fixed_points.FixedPointContinuation` class for more details
            about the available AUTO parameters.

        """

        logger.info('Starting the computation of the fixed points bifurcation diagram with model '+str(self.model_name))

        if self.fp_computed:
            logger.warning('Fixed point bifurcation diagram already computed. Aborting.')
            return None

        if initial_points is not None:
            self.initial_points = initial_points

        if 'NMX' not in continuation_kwargs:
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

            fp = FixedPointContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
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

            branch_order = 0

            parent_branch_number = sorted(self.fp_branches.keys())[branch_order]

            while True:

                if parent_branch_number not in self.fp_branches_with_all_bp_computed:

                    branch = self.fp_branches[parent_branch_number]
                    parent_continuation = branch['continuation']

                    logger.info('Continuing branching points of branch: ' + str(parent_branch_number))

                    forward_branching_points = parent_continuation.get_filtered_solutions_list(labels='BP', forward=True)
                    backward_branching_points = parent_continuation.get_filtered_solutions_list(labels='BP', forward=False)

                    direction_and_branching_points = ((1, forward_branching_points), (-1, backward_branching_points))

                    for direction, branching_points in direction_and_branching_points:

                        if direction == 1:
                            logger.info('Treating forward direction.')
                        else:
                            logger.info('Treating backward direction.')

                        for ibp, bp in enumerate(branching_points):

                            logger.info('Continuing fixed point out of branching point ' + str(ibp + 1))
                            logger.debug('First checking if branching point is not already computed in another branch.')
                            for bn in self.computed_bp_by_fp_branch:
                                found_solution = False
                                for ibpt in self.computed_bp_by_fp_branch[bn]:
                                    if ibpt >= 0:
                                        bpt = self.get_continuation(bn).get_solution_by_label('BP' + str(ibpt))
                                    else:
                                        bpt = self.get_continuation(bn).get_solution_by_label('-BP' + str(-ibpt))
                                    if self._check_if_solutions_are_close(bp, bpt, extra_comparison_parameters, comparison_tol):
                                        logger.debug('Branching point was already computed in branching point ' + str(ibpt) + ' of branch ' + str(bn) + '.')
                                        found_solution = True
                                        break
                                if found_solution:
                                    logger.debug('Skipping this point.')
                                    break
                            else:
                                logger.debug('Point is acceptable for continuation. Launching AUTO...')
                                used_continuation_kwargs = deepcopy(continuation_kwargs)
                                used_continuation_kwargs['ISW'] = -1
                                used_continuation_kwargs['PAR'] = {}
                                used_continuation_kwargs['IBR'] = br_num
                                fp = FixedPointContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
                                fp.make_continuation(bp, **used_continuation_kwargs)

                                logger.debug('Continuation done. Checking now against previous continuation...')
                                self._check_fp_continuation_against_itself(ncomp, fp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)
                                valid_branch = self._check_fp_continuation_against_other_fp_branches(ncomp, fp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol)

                                if valid_branch:
                                    self.fp_branches[abs(fp.branch_number)] = {'parameters': bp.PAR, 'continuation': fp, 'continuation_kwargs': used_continuation_kwargs}
                                    self.fp_parent[abs(fp.branch_number)] = parent_branch_number
                                    logger.info('Saving valid branch ' + str(br_num) + ' emanating from branch ' + str(parent_branch_number) + '.')
                                    br_num += 1
                                else:
                                    logger.info('Branching point ' + str(ibp + 1) + ' from branch ' + str(parent_branch_number) +
                                                ' resulted in non-valid branch. Not saving result.')
                                if parent_branch_number not in self.computed_bp_by_fp_branch:
                                    self.computed_bp_by_fp_branch[parent_branch_number] = list()
                                self.computed_bp_by_fp_branch[parent_branch_number].append(direction * (ibp+1))
                                ncomp += 1

                    self.fp_branches_with_all_bp_computed.append(parent_branch_number)
                    nrecomp += 1

                branch_order += 1
                try:
                    parent_branch_number = sorted(self.fp_branches.keys())[branch_order]
                except IndexError:
                    break

            if nrecomp == 0:
                break

        logger.info('Fixed points bifurcation diagram computation is over.')
        logger.info('All possible fixed point branches have been computed.')
        self.fp_computed = True

    def compute_periodic_orbits_diagram(self, end_level=10, extra_comparison_parameters=None, comparison_tol=2.e-2,
                                        remove_dubious_bp=True, max_number_bp=None, max_number_bp_detected=None, backward_bp_continuation=False, **continuation_kwargs):

        logger.info('Starting the computation of the periodic orbits bifurcation diagram with model '+str(self.model_name))
        logger.info('Computing periodic orbits up to level '+str(end_level))

        if not self.fp_computed:
            logger.warning('Fixed points diagram not computed. No initial data to start with.\n'
                           'Aborting...')
            return None

        if 'NMX' not in continuation_kwargs:
            continuation_kwargs['NMX'] = 9000
            logger.info('NMX parameters was not set, so setting it to 9000 points.')

        br_num = max(self.fp_branches.keys()) + 1
        self.level_reached = 0

        logger.info('First, beginning computation of the periodic orbits from Hopf bifurcations.')
        for parent_branch_number, branch in self.fp_branches.items():

            logger.info('Computing Hopf bifurcations of branch ' + str(parent_branch_number))

            forward_hb_list = branch['continuation'].get_filtered_solutions_list(labels='HB', forward=True)
            backward_hb_list = branch['continuation'].get_filtered_solutions_list(labels='HB', forward=False)

            direction_and_hopf_points = ((1, forward_hb_list), (-1, backward_hb_list))

            for direction, hb_list in direction_and_hopf_points:

                if direction == 1:
                    logger.info('Treating forward direction.')
                else:
                    logger.info('Treating backward direction.')

                for ihb, hb in enumerate(hb_list):

                    logger.info('Continuing PO out of Hopf point ' + str(ihb + 1) + '. Launching AUTO...')

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

                    hp = PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
                    hp.make_continuation(hb, max_bp=max_number_bp_detected, **used_continuation_kwargs)

                    logger.debug('Continuation done. Checking now against previous continuation...')
                    self._check_po_continuation_against_itself(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol, max_number_bp_detected)

                    self._check_po_continuation_against_other_fp_branches(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol, max_number_bp_detected)
                    valid_branch = self._check_po_continuation_against_other_po_branches(hp, used_continuation_kwargs, extra_comparison_parameters, comparison_tol, max_number_bp_detected)

                    if valid_branch:
                        self.po_branches[abs(hp.branch_number)] = {'parameters': hb.PAR, 'continuation': hp, 'continuation_kwargs': used_continuation_kwargs}
                        self.po_parent[abs(hp.branch_number)] = parent_branch_number
                        logger.info('Saving valid branch ' + str(br_num) + ' emanating from branch ' + str(parent_branch_number) + ' (Hopf point '+str(ihb+1)+').')
                        br_num += 1
                    else:
                        logger.info('Hopf point '+str(ihb+1)+' from branch '+str(parent_branch_number)+' resulted in non-valid branch. Not saving result.')

                    if parent_branch_number not in self.computed_hb_by_fp_branch:
                        self.computed_hb_by_fp_branch[parent_branch_number] = list()
                    self.computed_hb_by_fp_branch[parent_branch_number].append(direction * (ihb + 1))

            self.fp_branches_with_all_hb_computed.append(parent_branch_number)

        logger.info('Continuation of the periodic orbits from Hopf bifurcations has ended.')
        self.level_reached += 1
        if self.level_reached == end_level:
            logger.info('As demanded, finishing computation at level '+str(self.level_reached)+' ...')
            return

        logger.info('Beginning computation of the periodic orbits from detected branching and period doubling points.')
        self.restart_periodic_orbits_diagram(end_level=end_level, extra_comparison_parameters=extra_comparison_parameters,
                                             comparison_tol=comparison_tol,
                                             remove_dubious_bp=remove_dubious_bp, max_number_bp=max_number_bp,
                                             max_number_bp_detected=max_number_bp_detected,
                                             backward_bp_continuation=backward_bp_continuation, restart=False, **continuation_kwargs)

    def restart_periodic_orbits_diagram(self, end_level=10, extra_comparison_parameters=None, comparison_tol=2.e-2,
                                        remove_dubious_bp=True, max_number_bp=None, max_number_bp_detected=None,
                                        backward_bp_continuation=False, restart=True, **continuation_kwargs):

        if self.po_computed:
            logger.info('Computation up to level ' + str(end_level) + ' was asked but bifurcation diagram is complete.')
            logger.info('Nothing more to compute.')
        elif self.level_reached >= end_level:
            logger.info('Computation up to level ' + str(end_level) + ' was asked, but current level (' + str(self.level_reached) +
                        ') is already equal or above.')
            logger.info('Nothing more to compute.')
        elif self.level_reached == 0:
            logger.info('User asked for a restart, but actual bifurcation diagram level is 0 !')
            logger.info('Starting a new periodic orbit diagram...')
            self.compute_periodic_orbits_diagram(end_level=end_level, extra_comparison_parameters=extra_comparison_parameters,
                                                 remove_dubious_bp=remove_dubious_bp, max_number_bp=max_number_bp,
                                                 max_number_bp_detected=max_number_bp_detected,
                                                 backward_bp_continuation=backward_bp_continuation, **continuation_kwargs)
        else:
            if restart:
                logger.info('Restarting the computation of the periodic orbits from detected branching and period doubling points.')
                logger.info('Computing periodic orbits up to level ' + str(end_level))

            if 'NMX' not in continuation_kwargs:
                continuation_kwargs['NMX'] = 9000
                logger.info('NMX parameters was not set, so setting it to 9000 points.')

            br_num = max(self.po_branches.keys()) + 1

            while True:
                nrecomp = 0

                branch_order = 0
                max_branch_order_in_level = len(self.po_branches)
                parent_branch_number = sorted(self.po_branches.keys())[branch_order]

                logger.info('Entering level ' + str(self.level_reached + 1) + ' of continuation...')

                while branch_order <= max_branch_order_in_level:

                    if parent_branch_number not in self.po_branches_with_all_bp_computed:

                        branch = self.po_branches[parent_branch_number]
                        parent_continuation = branch['continuation']

                        logger.info('Continuing branching points of branch: ' + str(parent_branch_number))

                        forward_branching_points = parent_continuation.get_filtered_solutions_list(labels='BP', forward=True)
                        backward_branching_points = parent_continuation.get_filtered_solutions_list(labels='BP', forward=False)

                        if max_number_bp is not None:
                            forward_branching_points = forward_branching_points[:max_number_bp]
                            backward_branching_points = backward_branching_points[:max_number_bp]

                        direction_and_branching_points = ((1, forward_branching_points), (-1, backward_branching_points))

                        for direction, branching_points in direction_and_branching_points:

                            if direction == 1:
                                logger.info('Treating forward direction.')
                            else:
                                logger.info('Treating backward direction.')

                            for ibp, bp in enumerate(branching_points):
                                logger.info('Continuing PO of branching point '+str(ibp+1))
                                logger.debug('First checking if branching point is not already computed in another branch.')
                                for bn in self.computed_bp_by_po_branch:
                                    found_solution = False
                                    for ibpt in self.computed_bp_by_po_branch[bn]:
                                        if ibpt >= 0:
                                            bpt = self.get_continuation(bn).get_solution_by_label('BP' + str(ibpt))
                                        else:
                                            bpt = self.get_continuation(bn).get_solution_by_label('-BP' + str(-ibpt))
                                        if self._check_if_solutions_are_close(bp, bpt, extra_comparison_parameters, comparison_tol):
                                            logger.debug('Branching point was already computed in branching point ' + str(ibpt) + ' of branch '+str(bn)+'.')
                                            found_solution = True
                                            break
                                    if found_solution:
                                        logger.debug('Skipping this point.')
                                        break
                                else:
                                    logger.debug('Now checking the stability of the point.')
                                    try:
                                        bp_stability = np.array(parent_continuation.orbit_stability(direction * bp['PT']))
                                        # max_accept = 1. / np.nanmin(np.abs(np.where(bp_stability == 0, np.nan, bp_stability)))
                                        max_accept = np.finfo(np.float64).max * 1.e-10
                                        looks_dubious = np.max(bp_stability) > max_accept
                                    except ValueError:
                                        par_lst = parent_continuation.continuation_parameters
                                        par_val = [bp.PAR[p] for p in parent_continuation.continuation_parameters]
                                        ini_msg = str(par_lst) + " = " + str(par_val) + ' (branch ' + str(parent_branch_number) + ')'
                                        logger.error('No stability information found for PO point at ' + ini_msg + '. Something is wrong, not doing the continuation.')
                                        looks_dubious = True

                                    if looks_dubious and remove_dubious_bp:
                                        par_lst = parent_continuation.continuation_parameters
                                        par_val = [bp.PAR[p] for p in parent_continuation.continuation_parameters]
                                        ini_msg = str(par_lst) + " = " + str(par_val) + ' (branch ' + str(parent_branch_number) + ')'
                                        try:
                                            s = str(np.max(bp_stability))
                                        except ValueError:
                                            s = '[ unknown ]'
                                        logger.debug('Not saving results of PO point at ' + ini_msg + ' because it looks dubious. (max Floquet: ' + s + ' ).'
                                                     '\nSkipping to next one.')
                                        valid_branch = False
                                    else:
                                        logger.debug('Point is acceptable for continuation. Launching AUTO...')
                                        used_continuation_kwargs = deepcopy(self.po_branches[parent_branch_number]['continuation_kwargs'])
                                        used_continuation_kwargs['ISW'] = -1
                                        used_continuation_kwargs['NMX'] = continuation_kwargs['NMX']

                                        for param in continuation_kwargs:
                                            used_continuation_kwargs[param] = continuation_kwargs[param]

                                        if 'PAR' not in continuation_kwargs and 'PAR' in used_continuation_kwargs:
                                            _ = used_continuation_kwargs.pop('PAR')

                                        used_continuation_kwargs['IBR'] = br_num
                                        hp = PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
                                        hp.make_continuation(bp, only_forward=not backward_bp_continuation, max_bp=max_number_bp_detected, **used_continuation_kwargs)
                                        logger.debug('Continuation done. Checking now against previous continuation...')
                                        self._check_po_continuation_against_itself(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                                                   comparison_tol, max_number_bp_detected)

                                        self._check_po_continuation_against_other_fp_branches(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                                                              comparison_tol, max_number_bp_detected)
                                        valid_branch = self._check_po_continuation_against_other_po_branches(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                                                                             comparison_tol, max_number_bp_detected)

                                    if valid_branch:
                                        self.po_branches[abs(hp.branch_number)] = {'parameters': bp.PAR, 'continuation': hp, 'continuation_kwargs': used_continuation_kwargs}
                                        self.po_parent[abs(hp.branch_number)] = parent_branch_number
                                        if parent_branch_number not in self.valid_bp_by_po_branch:
                                            self.valid_bp_by_po_branch[parent_branch_number] = list()
                                        self.valid_bp_by_po_branch[parent_branch_number].append(direction * (ibp + 1))
                                        logger.info('Saving valid branch ' + str(br_num) + ' emanating from branch ' + str(parent_branch_number) +
                                                    ' (branching point '+str(ibp + 1)+').')
                                        br_num += 1
                                    else:
                                        logger.info('Branching point ' + str(ibp + 1) + ' from branch ' + str(parent_branch_number) +
                                                    ' resulted in non-valid branch. Not saving result.')
                                    if parent_branch_number not in self.computed_bp_by_po_branch:
                                        self.computed_bp_by_po_branch[parent_branch_number] = list()
                                    self.computed_bp_by_po_branch[parent_branch_number].append(direction * (ibp + 1))

                        self.po_branches_with_all_bp_computed.append(parent_branch_number)
                        nrecomp += 1

                        logger.info('Computation of branching points of branch: ' + str(parent_branch_number) + ' has ended.')

                    if parent_branch_number not in self.po_branches_with_all_pd_computed:

                        branch = self.po_branches[parent_branch_number]
                        parent_continuation = branch['continuation']

                        logger.info('Continuing period doubling points of branch: ' + str(parent_branch_number))

                        forward_period_doubling_points = parent_continuation.get_filtered_solutions_list(labels='PD', forward=True)
                        backward_period_doubling_points = parent_continuation.get_filtered_solutions_list(labels='PD', forward=False)

                        direction_and_period_doubling_points = ((1, forward_period_doubling_points), (-1, backward_period_doubling_points))

                        for direction, period_doubling_points in direction_and_period_doubling_points:

                            if direction == 1:
                                logger.info('Treating forward direction.')
                            else:
                                logger.info('Treating backward direction.')

                            for ipd, pd in enumerate(period_doubling_points):
                                logger.info('Continuing PO of period doubling point '+str(ipd+1))
                                logger.debug('First checking if period doubling point is not already computed in another branch.')

                                for bn in self.computed_pd_by_po_branch:
                                    found_solution = False
                                    for ipdt in self.computed_pd_by_po_branch[bn]:
                                        if ipdt >= 0:
                                            pdt = self.get_continuation(bn).get_solution_by_label('PD' + str(ipdt))
                                        else:
                                            pdt = self.get_continuation(bn).get_solution_by_label('-PD' + str(-ipdt))
                                        if self._check_if_solutions_are_close(pd, pdt, extra_comparison_parameters,
                                                                              comparison_tol):
                                            logger.debug('Period doubling point was already computed in period doubling point ' + str(ipdt) +
                                                         ' of branch '+str(bn)+'.')
                                            found_solution = True
                                            break
                                    if found_solution:
                                        logger.debug('Skipping this point.')
                                        break
                                else:
                                    logger.debug('Point is acceptable for continuation. Launching AUTO...')
                                    used_continuation_kwargs = deepcopy(self.po_branches[parent_branch_number]['continuation_kwargs'])
                                    used_continuation_kwargs['ISW'] = -1
                                    used_continuation_kwargs['NMX'] = continuation_kwargs['NMX']

                                    for param in continuation_kwargs:
                                        used_continuation_kwargs[param] = continuation_kwargs[param]

                                    if 'PAR' not in continuation_kwargs and 'PAR' in used_continuation_kwargs:
                                        _ = used_continuation_kwargs.pop('PAR')

                                    used_continuation_kwargs['IBR'] = br_num
                                    hp = PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
                                    hp.make_continuation(pd, max_bp=max_number_bp_detected, **used_continuation_kwargs)
                                    logger.debug('Continuation done. Checking now against previous continuation...')
                                    self._check_po_continuation_against_itself(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                                               comparison_tol, max_number_bp_detected)

                                    self._check_po_continuation_against_other_fp_branches(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                                                          comparison_tol, max_number_bp_detected)
                                    valid_branch = self._check_po_continuation_against_other_po_branches(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                                                                         comparison_tol, max_number_bp_detected)

                                    if valid_branch:
                                        self.po_branches[abs(hp.branch_number)] = {'parameters': pd.PAR, 'continuation': hp, 'continuation_kwargs': used_continuation_kwargs}
                                        self.po_parent[abs(hp.branch_number)] = parent_branch_number
                                        if parent_branch_number not in self.valid_pd_by_po_branch:
                                            self.valid_pd_by_po_branch[parent_branch_number] = list()
                                        self.valid_pd_by_po_branch[parent_branch_number].append(direction * (ipd + 1))
                                        logger.info('Saving branch ' + str(br_num) + ' emanating from branch ' + str(parent_branch_number) +
                                                    ' (period doubling point '+str(ipd + 1)+').')
                                        br_num += 1
                                    else:
                                        logger.info('Period doubling point ' + str(ipd + 1) + ' from branch ' + str(parent_branch_number) +
                                                    ' resulted in non-valid branch. Not saving result.')

                                    if parent_branch_number not in self.computed_pd_by_po_branch:
                                        self.computed_pd_by_po_branch[parent_branch_number] = list()
                                    self.computed_pd_by_po_branch[parent_branch_number].append(direction * (ipd + 1))

                        self.po_branches_with_all_pd_computed.append(parent_branch_number)
                        nrecomp += 1

                        logger.info('Computation of period doubling points of branch: ' + str(parent_branch_number) + ' has ended.')

                    branch_order += 1
                    try:
                        parent_branch_number = sorted(self.po_branches.keys())[branch_order]
                    except IndexError:
                        break

                if nrecomp == 0:
                    logger.info('No more solutions to continue, finishing ...')
                    break

                self.level_reached += 1
                logger.info('Computation of level ' + str(self.level_reached) + ' of continuation has ended.')
                if self.level_reached == end_level:
                    logger.info('As demanded, finishing computation at level ' + str(self.level_reached) + ' ...')
                    break

            logger.info('Finished computation at level ' + str(self.level_reached))
            for bn in self.po_branches:
                if bn not in self.po_branches_with_all_bp_computed or bn not in self.po_branches_with_all_pd_computed:
                    break
            else:
                self.po_computed = True
                logger.info('All possible periodic orbit branches have been computed.')

    def add_periodic_orbit(self, initial_data, extra_comparison_parameters=None, comparison_tol=2.e-2, max_number_bp_detected=None, only_forward=False,
                           **continuation_kwargs):

        logger.info('Continuing manually PO provided by user')

        if 'NMX' not in continuation_kwargs:
            continuation_kwargs['NMX'] = 9000
            logger.info('NMX parameters was not set, so setting it to 9000 points.')

        if self.po_branches:
            br_num = max(self.po_branches.keys()) + 1
        elif self.fp_branches:
            br_num = max(self.fp_branches.keys()) + 1
        else:
            br_num = 1

        logger.debug('First checking if branching point is not already computed in another branch.')
        found_solution = False
        # maybe tests on AUTOSolution must be put in try blocks (PO vs FP)
        if isinstance(initial_data, AUTOSolution):
            for bn in self.computed_bp_by_po_branch:
                for ibpt in self.computed_bp_by_po_branch[bn]:
                    if ibpt >= 0:
                        bpt = self.get_continuation(bn).get_solution_by_label('BP' + str(ibpt))
                    else:
                        bpt = self.get_continuation(bn).get_solution_by_label('-BP' + str(-ibpt))
                    if self._check_if_solutions_are_close(initial_data, bpt, extra_comparison_parameters, comparison_tol):
                        logger.debug('Branching point was already computed in branching point ' + str(ibpt) + ' of branch ' + str(bn) + '.')
                        found_solution = True
                        break
                if found_solution:
                    logger.debug('Skipping this point.')
                    break
            for bn in self.computed_pd_by_po_branch:
                found_solution = False
                for ipdt in self.computed_pd_by_po_branch[bn]:
                    if ipdt >= 0:
                        pdt = self.get_continuation(bn).get_solution_by_label('PD' + str(ipdt))
                    else:
                        pdt = self.get_continuation(bn).get_solution_by_label('-PD' + str(-ipdt))
                    if self._check_if_solutions_are_close(initial_data, pdt, extra_comparison_parameters,
                                                          comparison_tol):
                        logger.debug('Period doubling point was already computed in period doubling point ' + str(ipdt) +
                                     ' of branch ' + str(bn) + '.')
                        found_solution = True
                        break
                if found_solution:
                    logger.debug('Skipping this point.')
                    break
            for bn in self.computed_hb_by_fp_branch:
                found_solution = False
                for ihbt in self.computed_hb_by_fp_branch[bn]:
                    if ihbt >= 0:
                        hbt = self.get_continuation(bn).get_solution_by_label('HB' + str(ihbt))
                    else:
                        hbt = self.get_continuation(bn).get_solution_by_label('-HB' + str(-ihbt))
                    if self._check_if_solutions_are_close(initial_data, hbt, extra_comparison_parameters,
                                                          comparison_tol):
                        logger.debug('Hopf point was already computed in Hopf point ' + str(ihbt) +
                                     ' of branch ' + str(bn) + '.')
                        found_solution = True
                        break
                if found_solution:
                    logger.debug('Skipping this point.')
                    break
        else:
            for bn in list(self.fp_branches.keys()) + list(self.po_branches.keys()):
                for psol in self.get_continuation(bn).full_solutions_list:
                    try:
                        found_solution == (psol.initial_data == initial_data)
                    except:
                        pass
                    if found_solution:
                        logger.debug('Point was already computed in branch ' + str(bn) + '.')
                if found_solution:
                    logger.debug('Skipping this point.')
                    break

        if not found_solution:
            logger.debug('Point is acceptable for continuation. Launching AUTO...')
            used_continuation_kwargs = deepcopy(continuation_kwargs)
            used_continuation_kwargs['IBR'] = br_num
            used_continuation_kwargs['IPS'] = 2
            hp = PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
            hp.make_continuation(initial_data, only_forward=only_forward, max_bp=max_number_bp_detected, **used_continuation_kwargs)
            logger.debug('Continuation done. Checking now against previous continuation...')
            self._check_po_continuation_against_itself(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                       comparison_tol, max_number_bp_detected)

            self._check_po_continuation_against_other_fp_branches(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                                  comparison_tol, max_number_bp_detected)
            valid_branch = self._check_po_continuation_against_other_po_branches(hp, used_continuation_kwargs, extra_comparison_parameters,
                                                                                 comparison_tol, max_number_bp_detected)

            if valid_branch:
                if isinstance(initial_data, AUTOSolution):
                    parameters = initial_data.PAR
                elif 'PAR' in used_continuation_kwargs:
                    parameters = used_continuation_kwargs['PAR']
                else:
                    parameters = None
                self.po_branches[abs(hp.branch_number)] = {'parameters': parameters, 'continuation': hp, 'continuation_kwargs': used_continuation_kwargs}
                self.po_parent[abs(hp.branch_number)] = None
                logger.info('Saving valid branch ' + str(br_num) + ' emanating from user-provided data.')
                self.po_computed = False
                logger.info('There might be new periodic orbit branches to compute. '
                            'You might want to restart the automatic computation of the bifurcation diagram.')
            else:
                logger.info('User-provided data resulted in non-valid branch. Not saving result.')

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
        # store the pathname to pass to the updated class
        state['_path_name'] = self._path_name

        self.__dict__.clear()
        self.__dict__.update(state)

        filepath = 'c.' + self.model_name if self._path_name is None else os.path.join(self._path_name, 'c.' + self.model_name)
        self.config_object = ConfigParser(filepath)
        for branch_number in self.fp_branches:
            fp = FixedPointContinuation(self.model_name, self.config_object, path_name=self._path_name)
            fp.load('fp_' + str(branch_number) + '.pickle', load_initial_data=load_initial_data)
            self.fp_branches[branch_number]['continuation'] = fp
        for branch_number in self.po_branches:
            po = PeriodicOrbitContinuation(model_name=self.model_name, config_object=self.config_object, path_name=self._path_name)
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
        filepath = filename if self._path_name is None else os.path.join(self._path_name, filename)
        with open(filepath, 'wb') as f:
            pickle.dump(state, f, **kwargs)

    def load(self, filename, load_initial_data=True, **kwargs):
        try:
            filepath = filename if self._path_name is None else os.path.join(self._path_name, filename)
            with open(filepath, 'rb') as f:
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

            repeating, repeating_solutions = fp.check_for_repetitions(cpar_list, tol=tol, return_repeating_solutions=True, forward=True)
            recompute = False
            for i, v in enumerate(repeating[:-1]):
                if v and repeating[i + 1]:
                    recompute = True
                    break

            if recompute:

                first_repeating_sol = repeating_solutions[0]
                nmx = first_repeating_sol['PT'] + 1
                continuation_kwargs['NMX'] = nmx
                logger.info('Not storing full results of initial point ' + str(ncomp) + ' because it repeats itself (forward).'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                fp.make_forward_continuation(initial_data, **continuation_kwargs)

        if fp.continuation['backward'] is not None:

            repeating, repeating_solutions = fp.check_for_repetitions(cpar_list, tol=tol, return_repeating_solutions=True, forward=False)
            recompute = False
            for i, v in enumerate(repeating[:-1]):
                if v and repeating[i + 1]:
                    recompute = True
                    break

            if recompute:

                first_repeating_sol = repeating_solutions[0]
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

    def _check_po_continuation_against_other_fp_branches(self, continuation, continuation_kwargs, extra_comparison_parameters, tol, max_bp):

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
                hp.make_forward_continuation(initial_data, max_bp=max_bp, **continuation_kwargs)

            crossing, sol = hp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol, return_solutions=True, forward=False,
                                                     solutions_types=self._comparison_solutions_types)
            if crossing and not self._check_if_solutions_are_close(initial_data, sol, extra_comparison_parameters, tol):
                nmx = sol['PT'] + 1
                logger.info('Not storing full results of PO point at ' + ini_msg + ' because it connects to branch ' + str(n) + '.'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                continuation_kwargs['NMX'] = nmx
                hp.make_backward_continuation(initial_data, max_bp=max_bp, **continuation_kwargs)

    def _check_po_continuation_against_itself(self, continuation, continuation_kwargs, extra_comparison_parameters, tol, max_bp):

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

            repeating, repeating_solutions = hp.check_for_repetitions(cpar_list, tol=tol, return_repeating_solutions=True, forward=True)
            recompute = False
            for i, v in enumerate(repeating[:-1]):
                if v and repeating[i + 1]:
                    recompute = True
                    break

            if recompute:
                first_repeating_sol = repeating_solutions[0]
                nmx = first_repeating_sol['PT'] + 1
                logger.info('Not storing full results of PO point at ' + ini_msg + ' because it repeats itself (forward).'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                continuation_kwargs['NMX'] = nmx
                hp.make_forward_continuation(initial_data, max_bp=max_bp, **continuation_kwargs)

        if hp.continuation['backward'] is not None:

            repeating, repeating_solutions = hp.check_for_repetitions(cpar_list, tol=tol, return_repeating_solutions=True, forward=False)
            recompute = False
            for i, v in enumerate(repeating[:-1]):
                if v and repeating[i + 1]:
                    recompute = True
                    break

            if recompute:
                first_repeating_sol = repeating_solutions[0]
                nmx = first_repeating_sol['PT'] + 1
                logger.info('Not storing full results of PO point at ' + ini_msg + ' because it repeats itself (backward).'
                            '\nSaving only the relevant part. NMX set to ' + str(nmx))
                continuation_kwargs['NMX'] = nmx
                hp.make_backward_continuation(initial_data, max_bp=max_bp, **continuation_kwargs)

    def _check_po_continuation_against_other_po_branches(self, continuation, continuation_kwargs, extra_comparison_parameters, tol, max_bp):

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
                        hp.make_forward_continuation(initial_data, max_bp=max_bp, **continuation_kwargs)

                if not merge_forward and not remake_continuation:
                    cross_forward, sol = hp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol,
                                                                  return_solutions=True, forward=True, solutions_types=self._comparison_solutions_types)

                    if cross_forward:
                        nmx = sol['PT'] + 1
                        logger.info('Not storing full results of PO point at ' + ini_msg + ' because it connects forward to branch ' + str(n) + '.'
                                    '\nSaving only the relevant part. NMX set to ' + str(nmx))
                        continuation_kwargs['NMX'] = nmx
                        hp.make_forward_continuation(initial_data, max_bp=max_bp, **continuation_kwargs)

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
                        hp.make_backward_continuation(initial_data, max_bp=max_bp, **continuation_kwargs)

                if not merge_backward or not remake_continuation:
                    cross_backward, sol = hp.branch_possibly_cross(psol['continuation'], cpar_list, tol=tol,
                                                                   return_solutions=True, forward=False, solutions_types=self._comparison_solutions_types)

                    if cross_backward:
                        nmx = sol['PT'] + 1
                        logger.info('Not storing full results of  PO point at ' + ini_msg + ' because it connects backward to branch ' + str(n) + '.'
                                    '\nSaving only the relevant part. NMX set to ' + str(nmx))
                        continuation_kwargs['NMX'] = nmx
                        hp.make_backward_continuation(initial_data, max_bp=max_bp, **continuation_kwargs)

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
            self.fp_branches[b]['continuation'].plot_branch_parts(variables, ax=ax, **kwargs)
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
            self._figure_legend_handles = new_handles + self._figure_legend_handles
        elif len(self._figure_legend_handles) == self.number_of_fp_branches + self.number_of_po_branches:
            tmp_handles = self._figure_legend_handles[self.number_of_fp_branches:]
            self._figure_legend_handles = new_handles + tmp_handles
        else:
            self._figure_legend_handles = list()
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
            self.fp_branches[b]['continuation'].plot_branch_parts_3D(variables, ax=ax, **kwargs)
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
        elif len(self._figure_3d_legend_handles) == self.number_of_fp_branches + self.number_of_po_branches:
            tmp_handles = self._figure_3d_legend_handles[:self.number_of_fp_branches]
            self._figure_3d_legend_handles = new_handles + tmp_handles
        else:
            self._figure_3d_legend_handles = list()
            self._figure_3d_legend_handles.extend(new_handles)

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
            self.po_branches[b]['continuation'].plot_branch_parts(variables, ax=ax, **kwargs)
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
        elif len(self._figure_legend_handles) == self.number_of_fp_branches + self.number_of_po_branches:
            self._figure_legend_handles = self._figure_legend_handles[:self.number_of_fp_branches]
            self._figure_legend_handles.extend(new_handles)
        else:
            self._figure_legend_handles = list()
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
            self.po_branches[b]['continuation'].plot_branch_parts_3D(variables, ax=ax, **kwargs)
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
        elif len(self._figure_3d_legend_handles) == self.number_of_fp_branches + self.number_of_po_branches:
            self._figure_3d_legend_handles = self._figure_3d_legend_handles[:self.number_of_fp_branches]
            self._figure_3d_legend_handles.extend(new_handles)
        else:
            self._figure_3d_legend_handles = list()
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

    def plot_single_po_branch_and_solutions(
            self,
            branch_number,
            parameter=None,
            diagram_variables=(1,),
            solutions_variables=(0, 1),
            axes=None,
            figsize=(10, 16),
            solutions_tol=0.01,
            cmap=None,
            branch_diagram_kwargs=None,
            solution_diagram_kwargs=None,
            ):
        '''
            Plots a single branch, and all of the stored solutions on a single plot
        '''

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)
        
        if axes is None:
            fig, axes = plt.subplots(2, 1, figsize=figsize)

        if branch_diagram_kwargs is None:
            branch_diagram_kwargs = dict()

        if solution_diagram_kwargs is None:
            solution_diagram_kwargs = dict()

        if 'plot_kwargs' not in branch_diagram_kwargs:
            branch_diagram_kwargs['plot_kwargs'] = dict()

        if 'plot_kwargs' not in solution_diagram_kwargs:
            solution_diagram_kwargs['plot_kwargs'] = dict()

        if parameter is None:
            if self.po_branches:
                n = next(iter(self.po_branches))
                # Defaults to first parameter
                parameter = self.po_branches[n]['continuation_kwargs']['ICP'][0]
            else:
                return None

        if cmap is None:
            if 'cmap' in solution_diagram_kwargs['plot_kwargs']:
                cmap = solution_diagram_kwargs['plot_kwargs']['cmap']
            else:
                cmap = plt.get_cmap('Reds')
                solution_diagram_kwargs['plot_kwargs']['cmap'] = cmap
            branch_diagram_kwargs['plot_kwargs']['color'] = cmap[0.5]
        else:
            branch_diagram_kwargs['plot_kwargs']['color'] = cmap(0.5)
            solution_diagram_kwargs['plot_kwargs']['cmap'] = cmap
            
        self.po_branches[branch_number]['continuation'].plot_branch_parts(variables=(parameter, diagram_variables[0]), ax=axes[0], plot_sol_points=True, cmap=cmap, **branch_diagram_kwargs)

        # Plot scatter on branch of parameter values
        
        axes[0].set_title('Bifurcation diagram - Branch ' + str(branch_number))

        # Plot solution
        hp = self.po_branches[branch_number]['continuation']
        hp.plot_solutions(solutions_variables, ax=axes[1], parameter=parameter, value=None, color_solutions=True, tol=solutions_tol, **solution_diagram_kwargs)

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
