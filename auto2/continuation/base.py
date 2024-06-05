from abc import ABC, abstractmethod
import os
import sys
import warnings
import pickle

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
from auto.AUTOExceptions import AUTORuntimeError


# TODO: - Check what happens if parameters are integers
#       - Implement load and save method
#       - Change continutation attribute to a dict with "forward" and "backward" entries


class Continuation(ABC):

    def __init__(self, model_name, config_object):
        self.config_object = config_object
        self.model_name = model_name
        self.continuation = list()
        self.branch_number = None
        self.initial_data = None
        self.auto_filename_suffix = ""

        # plots default behaviours
        self._default_marker = None
        self._default_markersize = None
        self._default_linestyle = None
        self._default_linewidth = None

    @abstractmethod
    def make_continuation(self, initial_data, auto_suffix="", only_forward=False, **continuation_kwargs):
        pass

    @abstractmethod
    def make_forward_continuation(self, initial_data, auto_suffix="", **continuation_kwargs):
        pass

    @abstractmethod
    def make_backward_continuation(self, initial_data, auto_suffix="", **continuation_kwargs):
        pass

    def _get_dict(self):
        state = self.__dict__.copy()
        state['continuation'] = [None, None]
        if not isinstance(self.initial_data, np.ndarray) and self.initial_data is not None:
            state['initial_data'] = {key: self.initial_data[key] for key in ['BR', 'PT', 'TY name', 'TY number', 'Label']}
        return state

    @abstractmethod
    def _set_from_dict(self, state, load_initial_data=True):
        pass

    def auto_save(self, auto_suffix):
        if self.continuation:
            if self.continuation[0] is not None:
                ac.save(self.continuation[0], auto_suffix + '_forward')
            if self.continuation[1] is not None:
                ac.save(self.continuation[1], auto_suffix + '_backward')
        self.auto_filename_suffix = auto_suffix

    def auto_load(self, auto_suffix):
        self.continuation = list()
        try:
            r = ac.loadbd(auto_suffix + '_forward')
            self.continuation.append(r)
        except (FileNotFoundError, OSError, AUTORuntimeError):
            self.continuation.append(None)

        try:
            r = ac.loadbd(auto_suffix + '_backward')
            self.continuation.append(r)
        except (FileNotFoundError, OSError, AUTORuntimeError):
            self.continuation.append(None)

        if self.continuation[0] is None and self.continuation[1] is None:
            warnings.warn('Files not found. Unable to load data.')
        else:
            self.auto_filename_suffix = auto_suffix

    def save(self, filename=None, auto_filename_suffix=None, **kwargs):
        if auto_filename_suffix is None:
            warnings.warn('No AUTO filename suffix set. Using a default one.')
            if self.isfixedpoint:
                self.auto_filename_suffix = "fp_"+str(self.branch_number)
            else:
                self.auto_filename_suffix = "po_" + str(self.branch_number)
        else:
            self.auto_filename_suffix = auto_filename_suffix
        self.auto_save(self.auto_filename_suffix)
        if filename is None:
            warnings.warn('No pickle filename prefix provided. Using a default one.')
            if self.isfixedpoint:
                filename = "fp_"+str(self.branch_number)+'.pickle'
            else:
                filename = "po_"+str(self.branch_number)+'.pickle'
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

    @property
    def isfixedpoint(self):
        if self.config_object.parameters_dict[11] not in self.available_variables and 11 not in self.available_variables:
            return True
        else:
            return False

    @property
    def isperiodicorbit(self):
        return not self.isfixedpoint

    @property
    def available_variables(self):
        if self.continuation:
            if self.continuation[0] is not None:
                return self.continuation[0].data[0].keys()
            elif self.continuation[1] is not None:
                return self.continuation[1].data[0].keys()
            else:
                return None
        else:
            return None

    def diagnostics(self):
        if self.continuation:
            s = ''
            if self.continuation[0] is not None:
                s += 'Forward\n-----------------------\n'
                for t in self.continuation[0].data[0].diagnostics:
                    s += t['Text']
            if self.continuation[1] is not None:
                s += 'Backward\n-----------------------\n'
                for t in self.continuation[1].data[0].diagnostics:
                    s += t['Text']
            return s
        else:
            return None

    def print_diagnostics(self):
        if self.continuation:
            s = self.diagnostics()
            print(s)

    def point_diagnostic(self, idx):
        if self.continuation is not None:
            if isinstance(idx, str):
                if idx[0] == '-':
                    idx = self.find_solution_index(idx)
                    if idx is not None:
                        return self.continuation[1].data[0].diagnostics[idx]['Text']
                    else:
                        warnings.warn('No backward branch to show the diagnostic for.')
                        return None
                else:
                    idx = self.find_solution_index(idx)
                    if idx is not None:
                        return self.continuation[0].data[0].diagnostics[idx]
                    else:
                        warnings.warn('No forward branch to show the diagnostic for.')
                        return None

            if idx >= 0:
                if self.continuation[0] is not None:
                    return self.continuation[0].data[0].diagnostics[idx]
                else:
                    warnings.warn('No forward branch to show the diagnostic for.')
                    return None
            else:
                if self.continuation[1] is not None:
                    return self.continuation[1].data[0].diagnostics[-idx]['Text']
                else:
                    warnings.warn('No backward branch to show the diagnostic for.')
                    return None
        else:
            return None

    def find_solution_index(self, lab):
        if self.continuation:
            if lab[0] == "-":
                if self.solutions_label['backward'] is not None:
                    lab = lab[1:]
                    tgt = lab[:2]
                    sn = int(lab[2:])
                    idx = 0
                    count = 0
                    for la in self.solutions_label['backward']:
                        if la == tgt:
                            count += 1
                        if count >= sn:
                            break
                        idx += 1
                    else:
                        warnings.warn('No solution found.')
                        return None
                    return self.solutions_index['backward'][idx]
                else:
                    return None
            else:
                if self.solutions_label['forward'] is not None:
                    tgt = lab[:2]
                    sn = int(lab[2:])
                    idx = 0
                    count = 0
                    for la in self.solutions_label['forward']:
                        if la == tgt:
                            count += 1
                        if count >= sn:
                            break
                        idx += 1
                    else:
                        warnings.warn('No solution found.')
                        return None
                    return self.solutions_index['forward'][idx]
                else:
                    return None
        else:
            return None

    @property
    def stability(self):
        if self.continuation:
            s = list()
            if self.continuation[0] is not None:
                s.append(self.continuation[0].data[0].stability())
            else:
                s.append(list())
            if self.continuation[1] is not None:
                s.append(self.continuation[1].data[0].stability())
            else:
                s.append(list())
            return s
        else:
            return None

    @property
    def number_of_points(self):
        if self.continuation:
            n = list()
            if self.continuation[0] is not None:
                n.append(self.continuation[0].data[0].stability()[1])
            else:
                n.append(0)
            if self.continuation[1] is not None:
               n.append(self.continuation[1].data[0].stability()[1])
            else:
                n.append(0)
            return n
        else:
            return None

    @property
    def continuation_parameters(self):
        return self.continuation[0].c['ICP']

    @property
    def solutions_index(self):
        if self.continuation:
            d = dict()
            if self.continuation[0] is not None:
                idx = self.continuation[0].data[0].labels.getIndices()
                d['forward'] = idx
            else:
                d['forward'] = list()
            if self.continuation[1] is not None:
                idx = self.continuation[1].data[0].labels.getIndices()
                d['backward'] = idx
            else:
                d['backward'] = list()
            return d
        else:
            return None

    @property
    def solutions_label(self):
        indices = self.solutions_index
        if indices is not None:
            d = dict()
            if indices['forward']:
                idx = indices['forward']
                sl = list()
                for i in idx:
                    sl.append(str(self.continuation[0].data[0].labels.by_index[i].keys()).split("'")[1])
                d['forward'] = sl
            else:
                d['forward'] = list()
            if indices['backward']:
                idx = indices['backward']
                sl = list()
                for i in idx:
                    sl.append(str(self.continuation[1].data[0].labels.by_index[i].keys()).split("'")[1])
                d['backward'] = sl
            else:
                d['backward'] = list()
            return d
        else:
            return None

    @property
    def number_of_solutions(self):
        sl = list()
        if self.continuation:
            if self.continuation[0] is not None:
                sl.append(self.continuation[0].data[0].getLabels()[-1])
            else:
                sl.append(0)
            if self.continuation[1] is not None:
                sl.append(self.continuation[1].data[0].getLabels()[-1])
            else:
                sl.append(0)
            return sl
        else:
            return None

    @property
    def full_solutions_list_by_label(self):
        sd = self.solutions_list_by_direction_and_by_label
        if sd['backward']:
            sl = sd['backward']
        else:
            sl = list()
        if sd['forward']:
            sl.extend(sd['forward'])
        return sl

    @property
    def full_solutions_list(self):
        sd = self.solutions_list_by_direction
        if sd['backward']:
            if self.solutions_label['backward'][0] == 'EP':
                sl = sd['backward'][-1:0:-1]
            else:
                sl = sd['backward'][::-1]
        else:
            sl = list()
        if sd['forward']:
            sl.extend(sd['forward'])
        return sl

    @property
    def solutions_list_by_direction_and_by_label(self):
        sd = dict()
        if self.solutions_label is not None:
            sd['forward'] = list()
            sd['backward'] = list()
            if self.continuation[0] is not None:
                lab_list = list()
                for lab in self.solutions_label['forward']:
                    if lab not in lab_list:
                        lab_list.append(lab)
                for lab in lab_list:
                    sd['forward'].extend(self.continuation[0].getLabel(lab))
            if self.continuation[1] is not None:
                lab_list = list()
                for lab in self.solutions_label['backward']:
                    if lab not in lab_list:
                        lab_list.append(lab)
                for lab in lab_list:
                    sd['backward'].extend(self.continuation[1].getLabel(lab))
        return sd

    @property
    def solutions_list_by_direction(self):
        indices = self.solutions_index
        labels = self.solutions_label
        sd = dict()
        if indices is not None:
            sd['forward'] = list()
            sd['backward'] = list()
            if self.continuation[0] is not None:
                for lab, idx in zip(labels['forward'], indices['forward']):
                    if 'solution' in self.continuation[0].data[0].labels.by_index[idx][lab]:
                        sol = self.continuation[0].data[0].labels.by_index[idx][lab]['solution']
                        sd['forward'].append(sol)
            if self.continuation[1] is not None:
                for lab, idx in zip(labels['backward'], indices['backward']):
                    if 'solution' in self.continuation[1].data[0].labels.by_index[idx][lab]:
                        sol = self.continuation[1].data[0].labels.by_index[idx][lab]['solution']
                        sd['backward'].append(sol)
        return sd

    def solutions_parameters(self, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD', 'TR', 'LP'), forward=None):
        if not isinstance(parameters, (tuple, list)):
            parameters = [parameters]
        if isinstance(solutions_types, (tuple, list)):
            sl = self.get_filtered_solutions_list(labels=solutions_types, forward=forward)
        elif isinstance(solutions_types, str):
            if solutions_types == 'all':
                if forward is None:
                    sl = self.full_solutions_list
                elif forward:
                    sl = self.solutions_list_by_direction['forward']
                else:
                    sl = self.solutions_list_by_direction['backward']
        else:
            sl = self.full_solutions_list
        params = dict()
        for param in parameters:
            par_list = list()
            for s in sl:
                try:
                    par_list.append(max(s[param]))
                except TypeError:
                    par_list.append(float(s[param]))
            params[param] = par_list
        return np.squeeze(np.array(list(params.values()))).reshape((len(parameters), -1))

    def get_filtered_solutions_list(self, labels=None, indices=None, parameters=None, values=None, forward=None, tol=0.01):

        if parameters is not None:
            if isinstance(parameters, (list, tuple)):
                parameters = np.array(parameters)

            elif isinstance(parameters, str):
                parameters = np.array(parameters)[np.newaxis]

            if values is not None:

                if isinstance(values, (float, int)):
                    values = [values] * parameters.shape[0]
                    values = np.array(values)

                elif isinstance(values, (list, tuple)):
                    values = np.array(values)

                if len(values.shape) == 1:
                    values = values[:, np.newaxis]

                if values.shape[0] != parameters.shape[0]:
                    warnings.warn('Wrong number of values provided for the number of parameters test requested.')
                    return list()

            else:
                warnings.warn('No values provided for the parameters specified.')
                return list()

            if isinstance(tol, (list, tuple)):
                tol = np.array(tol)
            elif isinstance(tol, float):
                tol = np.array([tol] * values.shape[0])
            else:
                tol = np.array(tol)

        if forward is None:
            solutions_list = self.full_solutions_list
        elif forward:
            solutions_list = self.solutions_list_by_direction['forward']
        else:
            solutions_list = self.solutions_list_by_direction['backward']

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
        elif parameters is not None:
            new_solutions_list = list()
            for val in values.T:
                for sol in solutions_list:
                    sol_val = np.zeros_like(val)
                    for i, param in enumerate(parameters):
                        if isinstance(sol[param], np.ndarray):
                            sol_val[i] = max(sol[param])
                        else:
                            sol_val[i] = float(sol[param])
                    diff = sol_val - val
                    if np.all(np.abs(diff) < tol):
                        new_solutions_list.append(sol)
                        break
            solutions_list = new_solutions_list

        return solutions_list

    def plot_branche_parts(self, variables=(0, 1), ax=None, figsize=(10, 8), markersize=12., plot_kwargs=None, marker_kwargs=None,
                           excluded_labels=('UZ', 'EP'), variables_name=None):

        if not self.continuation:
            return None

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.gca()

        if plot_kwargs is None:
            plot_kwargs = dict()

        if marker_kwargs is None:
            marker_kwargs = dict()

        vars = self.available_variables

        for var in vars:
            try:
                if variables[0] in var:
                    var1 = var
                    break
            except:
                pass
        else:
            try:
                var1 = vars[variables[0]]
            except:
                var1 = vars[0]

        for var in vars:
            try:
                if variables[1] in var:
                    var2 = var
                    break
            except:
                pass
        else:
            try:
                var2 = vars[variables[1]]
            except:
                var2 = vars[1]

        for branche_part in [0, 1]:

            if self.continuation[branche_part] is not None:

                labels = list()
                for j, coords in enumerate(zip(self.continuation[branche_part][var1], self.continuation[branche_part][var2])):
                    lab = self.continuation[branche_part].data[0].getIndex(j)['TY name']
                    if lab == 'No Label':
                        pass
                    else:
                        labels.append((coords, lab))

                pidx = 0
                for idx in self.stability[branche_part]:
                    if idx < 0:
                        ls = '-'
                    else:
                        ls = '--'
                    plot_kwargs['linestyle'] = ls
                    lines_list = ax.plot(self.continuation[branche_part][var1][pidx:abs(idx)], self.continuation[branche_part][var2][pidx:abs(idx)], **plot_kwargs)
                    c = lines_list[0].get_color()
                    plot_kwargs['color'] = c
                    pidx = abs(idx)
                if excluded_labels != 'all':
                    for label in labels:
                        coords = label[0]
                        lab = label[1]
                        if lab not in excluded_labels:
                            ax.text(coords[0], coords[1], r'${\bf ' + lab + r'}$', fontdict={'family': 'sans-serif', 'size': markersize}, va='center', ha='center', **marker_kwargs,
                                    clip_on=True)

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

    def plot_branche_parts_3D(self, variables=(0, 1, 2), ax=None, figsize=(10, 8), markersize=12., plot_kwargs=None, marker_kwargs=None,
                           excluded_labels=('UZ', 'EP'), variables_name=None):

        if not self.continuation:
            return None

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(projection='3d')

        if plot_kwargs is None:
            plot_kwargs = dict()

        if marker_kwargs is None:
            marker_kwargs = dict()

        vars = self.available_variables

        for var in vars:
            try:
                if variables[0] in var:
                    var1 = var
                    break
            except:
                pass
        else:
            try:
                var1 = vars[variables[0]]
            except:
                var1 = vars[0]

        for var in vars:
            try:
                if variables[1] in var:
                    var2 = var
                    break
            except:
                pass
        else:
            try:
                var2 = vars[variables[1]]
            except:
                var2 = vars[1]

        for var in vars:
            try:
                if variables[2] in var:
                    var3 = var
                    break
            except:
                pass
        else:
            try:
                var3 = vars[variables[2]]
            except:
                var3 = vars[2]

        for branche_part in [0, 1]:

            if self.continuation[branche_part] is not None:

                labels = list()
                for j, coords in enumerate(zip(self.continuation[branche_part][var1], self.continuation[branche_part][var2], self.continuation[branche_part][var3])):
                    lab = self.continuation[branche_part].data[0].getIndex(j)['TY name']
                    if lab == 'No Label':
                        pass
                    else:
                        labels.append((coords, lab))

                pidx = 0
                for idx in self.stability[branche_part]:
                    if idx < 0:
                        ls = '-'
                    else:
                        ls = '--'
                    plot_kwargs['linestyle'] = ls
                    lines_list = ax.plot(self.continuation[branche_part][var1][pidx:abs(idx)], self.continuation[branche_part][var2][pidx:abs(idx)],
                                         self.continuation[branche_part][var3][pidx:abs(idx)], **plot_kwargs)
                    c = lines_list[0].get_color()
                    plot_kwargs['color'] = c
                    pidx = abs(idx)
                if excluded_labels != 'all':
                    for label in labels:
                        coords = label[0]
                        lab = label[1]
                        if lab not in excluded_labels:
                            ax.text(coords[0], coords[1], coords[2], r'${\bf ' + lab + r'}$', fontdict={'family': 'sans-serif', 'size': markersize}, va='center', ha='center',
                                    **marker_kwargs, clip_on=True)

        if variables_name is None:
            ax.set_xlabel(var1)
            ax.set_ylabel(var2)
            ax.set_zlabel(var3)
        else:
            if isinstance(variables_name, dict):
                ax.set_xlabel(variables_name[var1])
                ax.set_ylabel(variables_name[var2])
                ax.set_zlabel(variables_name[var3])
            else:
                ax.set_xlabel(variables_name[0])
                ax.set_ylabel(variables_name[1])
                ax.set_zlabel(variables_name[2])
        return ax

    def plot_solutions(self, variables=(0, 1), ax=None, figsize=(10, 8), markersize=None, marker=None, linestyle=None,
                       linewidth=None, plot_kwargs=None, labels=None, indices=None, parameter=None, value=None,
                       variables_name=None, tol=0.01):

        if markersize is None:
            markersize = self._default_markersize
        if marker is None:
            marker = self._default_marker
        if linestyle is None:
            linestyle = self._default_linestyle
        if linewidth is None:
            linewidth = self._default_linewidth

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.gca()

        if plot_kwargs is None:
            plot_kwargs = dict()

        solutions_list = self.get_filtered_solutions_list(labels, indices, parameter, value, tol)

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

    def plot_solutions_3D(self, variables=(0, 1, 2), ax=None, figsize=(10, 8), markersize=None, marker=None, linestyle=None,
                       linewidth=None, plot_kwargs=None, labels=None, indices=None, parameter=None, value=None,
                       variables_name=None, tol=0.01):

        if markersize is None:
            markersize = self._default_markersize
        if marker is None:
            marker = self._default_marker
        if linestyle is None:
            linestyle = self._default_linestyle
        if linewidth is None:
            linewidth = self._default_linewidth

        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(projection='3d')

        if plot_kwargs is None:
            plot_kwargs = dict()

        solutions_list = self.get_filtered_solutions_list(labels, indices, parameter, value, tol)

        vars = self.config_object.variables

        if variables[0] in vars:
            var1 = variables[0]
        else:
            try:
                var1 = vars[variables[0]]
            except:
                var1 = vars[0]

        if variables[1] in vars:
            var2 = variables[1]
        else:
            try:
                var2 = vars[variables[1]]
            except:
                var2 = vars[1]

        if variables[2] in vars:
            var3 = variables[2]
        else:
            try:
                var3 = vars[variables[2]]
            except:
                var3 = vars[2]

        for sol in solutions_list:
            x = sol[var1]
            y = sol[var2]
            z = sol[var3]
            line_list = ax.plot(x, y, z, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth,
                                **plot_kwargs)
            c = line_list[0].get_color()
            plot_kwargs['color'] = c

        if variables_name is None:
            ax.set_xlabel(var1)
            ax.set_ylabel(var2)
            ax.set_zlabel(var3)
        else:
            if isinstance(variables_name, dict):
                ax.set_xlabel(variables_name[var1])
                ax.set_ylabel(variables_name[var2])
                ax.set_zlabel(variables_name[var3])
            else:
                ax.set_xlabel(variables_name[0])
                ax.set_ylabel(variables_name[1])
                ax.set_zlabel(variables_name[2])
        return ax

    def same_solutions_as(self, other, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD', 'TR', 'LP'), tol=2.e-2):

        ssol = self.solutions_parameters(parameters, solutions_types)
        osol = other.solutions_parameters(parameters, solutions_types)

        npar = ssol.shape[0]
        if isinstance(tol, float):
            tol = [tol] * npar

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        ssol = _sort_arrays(ssol, npar, tol)
        osol = _sort_arrays(osol, npar, tol)

        if ssol.shape != osol.shape:
            return False
        else:
            diff = ssol - osol
            return np.all(np.abs(diff).T < tol)

    def solutions_in(self, other, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD', 'TR', 'LP'), tol=2.e-2, return_parameters=False, return_solutions=False):
        res, params, sol = self.solutions_part_of(other, parameters, solutions_types, tol, True, True, None)
        if res:
            res = [params.shape[1] == self.number_of_solutions]
        else:
            res = [res]
        if return_parameters:
            res.append(params)
        if return_solutions:
            res.append(sol)

        if len(res) == 1:
            return res[0]
        else:
            return res

    def solutions_part_of(self, other, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD', 'TR', 'LP'), tol=2.e-2, return_parameters=False, return_solutions=False, forward=None):

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        ssol = self.solutions_parameters(parameters, solutions_types, forward=forward)
        osol = other.solutions_parameters(parameters, solutions_types)

        npar = ssol.shape[0]
        if isinstance(tol, float):
            tol = [tol] * npar

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        idx_list = list()
        for i in range(ssol.shape[1]):
            for j in range(osol.shape[1]):
                dif = ssol[:, i] - osol[:, j]
                if np.all(np.abs(dif) < tol):
                    idx_list.append(i)
                    break

        res = [len(idx_list) > 1]
        if return_parameters:
            if res[0]:
                res.append(ssol[:, idx_list])
            else:
                res.append(np.empty(0))
        if return_solutions:
            if res[0]:
                res.append(self.get_filtered_solutions_list(parameters=parameters, values=ssol[:, idx_list], tol=tol, forward=forward))
            else:
                res.append(None)

        if len(res) == 1:
            return res[0]
        else:
            return res

    def branch_possibly_cross(self, other, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD', 'TR', 'LP'), tol=2.e-2, return_parameters=False, return_solutions=False, forward=None):

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        ssol = self.solutions_parameters(parameters, solutions_types, forward=forward)
        osol = other.solutions_parameters(parameters, solutions_types)

        npar = ssol.shape[0]
        if isinstance(tol, float):
            tol = [tol] * npar

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        idx_list = list()
        for i in range(ssol.shape[1]):
            for j in range(osol.shape[1]):
                dif = ssol[:, i] - osol[:, j]
                if np.all(np.abs(dif) < tol):
                    idx_list.append(i)
                    break

        res = [len(idx_list) == 1]
        if return_parameters:
            if res[0]:
                res.append(ssol[:, idx_list])
            else:
                res.append(np.empty(0))
        if return_solutions:
            if res[0]:
                res.append(self.get_filtered_solutions_list(parameters=parameters, values=ssol[:, idx_list], tol=tol, forward=forward)[0])
            else:
                res.append(None)

        if len(res) == 1:
            return res[0]
        else:
            return res

    def summary(self):

        summary_str = ""
        if self.continuation[0] is not None:
            summary_str += "Forward\n"
            summary_str += "=======\n"
            summary_str += self.continuation[0].summary() + "\n\n"

        if self.continuation[1] is not None:
            summary_str += "Backward\n"
            summary_str += "========\n"
            summary_str += self.continuation[1].summary()

        return summary_str

    def print_summary(self):

        if self.continuation:
            s = self.summary()
            print(s)

    def check_for_repetitions(self, parameters, tol=2.e-2, return_parameters=False, return_non_repeating_solutions=False,
                              return_repeating_solutions=False, forward=None):

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        valid_solution = dict()
        valid_solution['backward'] = list()
        valid_solution['forward'] = list()
        if forward:
            valid_solution['backward'] = None
        elif not forward:
            valid_solution['forward'] = None

        forward_dict = {'forward': True, 'backward': False}

        sol = dict()
        idx_dict = dict()
        ridx_dict = dict()
        for direction in ['forward', 'backward']:
            if valid_solution[direction] is not None:
                idx_dict[direction] = list()
                ridx_dict[direction] = list()
                idx_list = list()
                ssol = self.solutions_parameters(parameters, solutions_types='all', forward=forward_dict[direction])
                sol[direction] = ssol

                npar = ssol.shape[0]
                if isinstance(tol, float):
                    tol = [tol] * npar

                if isinstance(tol, (list, tuple)):
                    tol = np.array(tol)

                for i in range(ssol.shape[1]-1, -1, -1):
                    for j in range(i-1, 0, -1):
                        dif = ssol[:, i] - ssol[:, j]
                        if np.all(np.abs(dif) < tol):
                            idx_list.append(i)
                            break

                for i in range(ssol.shape[1]):
                    if i not in idx_list:
                        valid_solution[direction].append(True)
                        idx_dict[direction].append(i)

                    else:
                        valid_solution[direction].append(False)
                        ridx_dict[direction].append(i)

        if forward is None:
            if valid_solution['backward'] is not None:
                res = [valid_solution['backward']]
            else:
                res = list([])
            if valid_solution['forward'] is not None:
                res[0] += valid_solution['forward']
        elif forward:
            if valid_solution['forward'] is not None:
                res = [valid_solution['forward']]
            else:
                res = list([])
        else:
            if valid_solution['backward'] is not None:
                res = [valid_solution['backward']]
            else:
                res = list([])

        if return_parameters:
            if forward is None:
                params = np.concatenate((sol['backward'][:, idx_dict['backward']], sol['foward'][:, idx_dict['forward']]), axis=-1)
            elif forward:
                params = sol['forward'][:, idx_dict['forward']]
            else:
                params = sol['backward'][:, idx_dict['backward']]

            res.append(params)

        if return_non_repeating_solutions:
            if forward is None:
                sols = (self.get_filtered_solutions_list(parameters=parameters, values=sol['backward'][:, idx_dict['backward']], tol=tol) +
                        self.get_filtered_solutions_list(parameters=parameters, values=sol['forward'][:, idx_dict['forward']], tol=tol))
            elif forward:
                sols = self.get_filtered_solutions_list(parameters=parameters, values=sol['forward'][:, idx_dict['forward']], tol=tol)
            else:
                sols = self.get_filtered_solutions_list(parameters=parameters, values=sol['backward'][:, idx_dict['backward']], tol=tol)

            res.append(sols)

        if return_repeating_solutions:
            if forward is None:
                rsols = (self.get_filtered_solutions_list(parameters=parameters, values=sol['backward'][:, ridx_dict['backward']], tol=tol) +
                         self.get_filtered_solutions_list(parameters=parameters, values=sol['forward'][:, ridx_dict['forward']], tol=tol))
            elif forward:
                rsols = self.get_filtered_solutions_list(parameters=parameters, values=sol['forward'][:, ridx_dict['forward']], tol=tol)
            else:
                rsols = self.get_filtered_solutions_list(parameters=parameters, values=sol['backward'][:, ridx_dict['backward']], tol=tol)

            res.append(rsols)

        if len(res) == 1:
            return res[0]
        else:
            return res


def _sort_arrays(sol, npar, tol):

    srt = np.squeeze(np.argsort(np.ascontiguousarray(sol.T).view(','.join(['f8'] * npar)), order=['f' + str(i) for i in range(npar)], axis=0).T)

    ssol = sol[:, srt].reshape((npar, -1))

    for n in range(1, npar):
        while True:
            nc = 0
            for i in range(ssol.shape[1] - 1):
                if abs(ssol[n - 1, i + 1] - ssol[n - 1, i]) < tol[n - 1] and ssol[n, i + 1] < ssol[n, i]:
                    nc += 1
                    a = ssol[:, i + 1].copy()
                    b = ssol[:, i].copy()
                    ssol[:, i] = a
                    ssol[:, i+1] = b
            if nc == 0:
                break

    return ssol

