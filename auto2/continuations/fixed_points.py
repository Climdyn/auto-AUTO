"""

    Fixed Point Continuation class definition
    =========================================

    This module implements the fixed point continuation methods.

"""

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
    """Class for the fixed point continuations in auto-AUTO.

    Parameters
    ----------
    model_name: str
        The name of the model to load. Used to load the .f90 file of the provided model.
    config_object: ~auto2.parsers.config.ConfigParser
        A loaded ConfigParser object.
    path_name: str, optional
        The directory path where files are read/saved.
        If `None`, defaults to current working directory.

    Attributes
    ----------
    model_name: str
        The name of the loaded model. Used to load the .f90 file of the provided model.
    config_object: ~auto2.parsers.config.ConfigParser
        A loaded ConfigParser object.
    continuation: dict
        Dictionary holding the forward and backward continuation data.
    branch_number: int
        The AUTO branch number attributed to the continuation(s). This number can be set manually by passing the IBR AUTO parameter when
        starting the continuations (see the :meth:`make_continuation` documentation for more details).
    initial_data: ~numpy.ndarray or AUTO solution object
        The initial data used to start the continuation(s).
    auto_filename_suffix: str
        Suffix for the |AUTO| files used to save the continuation(s) data and parameters on disk.

    """

    def __init__(self, model_name, config_object, path_name=None):

        Continuation.__init__(self, model_name, config_object, path_name)

        # plots default behaviours
        self._default_marker = 'x'
        self._default_markersize = 6.
        self._default_linestyle = ' '
        self._default_linewidth = 1.2

    def make_continuation(self, initial_data, auto_suffix="", only_forward=False, **continuation_kwargs):
        """Make both forward and backward (if possible) continuation of fixed points.

        Parameters
        ----------
        initial_data: ~numpy.ndarray or AUTOSolution
            Initial data used to start the continuation(s).
        auto_suffix: str, optional
            Suffix to use for the |AUTO| and Pickle files used to save the continuation(s)
            data, parameters and metadata on disk. If not provided, does not save the data on disk.
        only_forward: bool, optional
            If `True`, compute only the forward continuation (positive `DS` parameter).
            If `False`, compute in both backward and forward direction.
            Default to `False`.
        continuation_kwargs: dict
            Keyword arguments to be passed to the |AUTO| continuation.
            See below for further details.

        Notes
        -----

        **AUTO Continuation Keyword Arguments**: See Section 10.8 in the |AUTO| documentation for further details.
        We provide below the most important ones:

        Other Parameters
        ----------------

        DS: float, optional
            AUTO uses pseudo-arclength continuation for following solution families.
            The pseudo-arclength stepsize is the distance between the current solution and the next solution on a family.
            By default, this distance includes all state variables (or state functions) and all free parameters.
            The constant `DS` defines the pseudo-arclength stepsize to be used for the first attempted step along any family.
            DS may be chosen positive or negative; changing its sign reverses the direction of computation.
            The relation `DSMIN` ≤ | `DS` | ≤ `DSMAX` must be satisfied. The precise choice of `DS` is problem-dependent.
        DSMIN: float, optional
            This is minimum allowable absolute value of the pseudo-arclength stepsize.
            `DSMIN` must be positive. It is only effective if the pseudo-arclength step is adaptive, i.e., if `IADS`>0.
            The choice of `DSMIN` is highly problem-dependent.
        DSMAX: float, optional
            The maximum allowable absolute value of the pseudo-arclength stepsize.
            `DSMAX` must be positive.
            It is only effective if the pseudo-arclength step is adaptive, i.e., if `IADS`>0.
            The choice of `DSMAX` is highly problem-dependent.
        NMX: int, optional
            The maximum number of steps to be taken along the branch.
        IBR: int, optional
            This constant specifies the initial branch number BR that is used. The default `IBR=0` means that
            that this number is determined automatically.
        ILP: int, optional
            * If `ILP=0`: No detection of folds. This is the recommended choice in the AUTO documentation.

            * If `ILP=1`: Detection of folds. To be used if subsequent fold continuation is intended.
        SP: list(str), optional
            This constant controls the detection of bifurcations and adds stopping conditions.
            It is specified as a list of bifurcation type strings followed by an optional number.
            If this number is `0`, then the detection of this bifurcation is turned off, and if it is missing then the detection is turned on.
            A number `n` greater than zero specifies that the continuation should stop as soon as the nth bifurcation of this type has been
            reached.
            Examples:

            - `SP=[’LP0’]` turn off detection of folds.

            - `SP=[’LP’,’HB3’,’BP0’,’UZ3’]` turn on the detection of folds and Hopf bifurcations, turn off detection of branch points
              and stop at the third Hopf bifurcation or third user defined point, whichever comes first.
        ISP: int, optional
            This constant controls the detection of Hopf bifurcations, branch points, period-doubling bifurcations, and torus bifurcations:

            * If `ISP=0` This setting disables the detection of Hopf bifurcations, branch points, period doubling bifurcations,
              and torus bifurcations and the computation of Floquet multipliers.

            * If `ISP=1` Branch points and Hopf bifurcations are detected for algebraic equations. Branch points, period-doubling
              bifurcations and torus bifurcations are not detected for periodic solutions and boundary value problems.
              However, Floquet multipliers are computed.

            * If `ISP=2` This setting enables the detection of all special solutions.
              For periodic solutions and rotations, the choice `ISP=2` should be used with care,
              due to potential inaccuracy in the computation of the linearized Poincar ́e map and possible rapid variation of the
              Floquet multipliers.

            * If `ISP=3` Hopf bifurcations will not be detected. Branch points will be detected, and AUTO will monitor
              the Floquet multipliers. Period-doubling and torus bifurcations will go undetected. This option is useful for
              certain problems with non-generic Floquet behavior.

            * If `ISP=4` Branch points and Hopf bifurcations are detected for algebraic equations. Branch points are not detected
              for periodic solutions and boundary value problems.
              AUTO will monitor the Floquet multipliers, and period-doubling and torus bifurcations will be detected.
        ISW: int, optional
            This constant controls branch switching at branch points for the case of differential equations. Note that branch switching
            is automatic for algebraic equations.

            * If `ISW=1` This is the normal value of `ISW`.

            * If `ISW=-1` If `IRS` is the label of a branch point or a period-doubling bifurcation then branch switching will be done.
              For period doubling bifurcations it is recommended that `NTST` be increased.

            * If `ISW=2` If IRS is the label of a fold, a Hopf bifurcation point, a period-doubling, a torus bifurcation,
              or, in a non-generic (symmetric) system, a branch point then a locus of such points will be computed.

            * If `ISW=3` If `IRS` is the label of a branch point in a generic (non-symmetric) system then a locus of such points will
              be computed. Two additional free parameters must be specified for such continuations
        MXBF: int, optional
            This constant, which is effective for algebraic problems only, sets the maximum number of bifurcations to be treated.
            Additional branch points will be noted, but the corresponding bifurcating families will not be computed.
        dat: str, optional
            This constant, where `dat=’filename’`, sets the name of a user-supplied ASCII data file `filename.dat`,
            from which the continuation is to be restarted.
            The first column in the data file denotes the time, which does not need to be rescaled to the
            interval `[0, 1]`, and further columns the coordinates of the solution. The parameter `IRS` must be
            set to `0`.
        PAR: dict, optional
             Determines the parameter values for the continuation to start from.
             Should be entred as `{'parameter_name': parameter_value}`.
        IRS: int or str, optional
            This constant sets the label of the solution where the computation is to be restarted.
            Setting `IRS=0` is typically used in the first run of a new problem.
            To restart the computation at the `n`-th label, use `IRS=n`.
            To restart the computation at a specific label (for example `HB12`), use `IRS=HB12`.
        TY: str, optional
            This constant modifies the type from the restart solution.
            This is sometimes useful in conservative or extended systems, declaring a regular point to be
            a Hopf bifurcation point `TY=’HB’` or a branch point `TY=’BP’`.
        IPS: int, optional
            This constant defines the problem type:

            * If `IPS=0` algebraic bifurcation problem.
            * If `IPS=1` stationary solutions of ODEs with detection of Hopf bifurcations.
            * If `IPS=-1` fixed points of the discrete dynamical systems.
            * If `IPS=-2` time integration using implicit Euler.
            * If `IPS=2`  computation of periodic solutions.
            * If `IPS=4`  boundary value problems.
            * If `IPS=5`  algebraic optimization problems.
            * If `IPS=7`  boundary value problem with computation of Floquet multipliers.
            * If `IPS=9`  option is used in connection with the HomCont algorithms
            * If `IPS=11` spatially uniform solutions of a system of parabolic PDEs, with detection of traveling wave bifurcations.
            * If `IPS=12` continuation of traveling wave solutions to a system of parabolic PDEs.
            * If `IPS=14` time evolution for a system of parabolic PDEs subject to periodic boundary conditions.

        """

        self.make_forward_continuation(initial_data, "", **continuation_kwargs)
        if not only_forward:
            self.make_backward_continuation(initial_data, "", **continuation_kwargs)

        if auto_suffix:
            self.auto_save(auto_suffix)

    def make_forward_continuation(self, initial_data, auto_suffix="", **continuation_kwargs):
        """Make the forward continuation of fixed points.

        Parameters
        ----------
        initial_data: ~numpy.ndarray or AUTOSolution
            Initial data used to start the continuation(s).
        auto_suffix: str, optional
            Suffix to use for the |AUTO| and Pickle files used to save the continuation(s)
            data, parameters and metadata on disk. If not provided, does not save the data on disk.
        continuation_kwargs: dict
            Keyword arguments to be passed to the |AUTO| continuation.
            See below for further details

        Notes
        -----

        **AUTO Continuation Keyword Arguments**: See Section 10.8 in the |AUTO| documentation for further details.
        The most important ones are provided in the documentation of the :meth:`make_continuation` method.

        """
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
        """Make the backward continuation of fixed points.

        Parameters
        ----------
        initial_data: ~numpy.ndarray or AUTOSolution
            Initial data used to start the continuation(s).
        auto_suffix: str, optional
            Suffix to use for the |AUTO| and Pickle files used to save the continuation(s)
            data, parameters and metadata on disk. If not provided, does not save the data on disk.
        continuation_kwargs: dict
            Keyword arguments to be passed to the |AUTO| continuation.
            See below for further details

        Notes
        -----

        **AUTO Continuation Keyword Arguments**: See Section 10.8 in the |AUTO| documentation for further details.
        The most important ones are provided in the documentation of the :meth:`make_continuation` method.

        """
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
        """Given a point index, returns the stability values (Eiganvalues) for this point.

        Parameters
        ----------
        idx: str, int
            |AUTO| index of the point to give the diagnostic from.
            If `idx`:

            * is a string, it should identify the label of the point, assuming the
              point has a label. Backward continuation label must be preceded by
              a `'-'` character. E.g. `'HB2'` would identify the second Hopf bifurcation
              point of the forward continuation, while `'-UZ3'` identifies the third
              user-defined label of the backward branch.
            * is an integer, it corresponds to the index of the point. Positive
              integers identify points on the forward continuation, while negative integers
              identify points on the backward continuation. The zero index identifies the
              initial point of the continuations.
            

        Returns
        -------
        ~numpy.ndarray
        """
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
        """bool: Whether this class stores fixed point continuations."""
        return True

    @property
    def isperiodicorbit(self):
        """bool: Whether this class stores periodic orbit continuations."""
        return False
