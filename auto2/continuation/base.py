import re
from abc import ABC, abstractmethod
import os
import sys
import warnings
import logging
import pickle
from contextlib import contextmanager

import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger('logger')

try:
    auto_directory = os.environ['AUTO_DIR']

    for path in sys.path:
        if auto_directory in path:
            break
    else:
        # sys.path.append(auto_directory + '/python/auto')
        sys.path.append(sys.path.join(auto_directory, 'python'))
except KeyError:
    logger.warning('Unable to find auto directory environment variable.')

import auto.AUTOCommands as ac
from auto.AUTOExceptions import AUTORuntimeError


# TODO: - Check what happens if parameters are integers


class Continuation(ABC):

    def __init__(self, model_name, config_object, path_name=None):
        self.config_object = config_object
        self.model_name = model_name
        self.continuation = dict()
        self.branch_number = None
        self.initial_data = None
        self.auto_filename_suffix = ""

        if path_name is None:
            self._path_name = None
        else:
            if os.path.exists(path_name):
                self._path_name = path_name
            else:
                warnings.warn("Path name given does not exist.")
                self._path_name = None

        # plots default behaviours
        self._default_marker = None
        self._default_markersize = None
        self._default_linestyle = None
        self._default_linewidth = None

        # options
        self._retry = 3

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
        state['continuation'] = {
            'forward': None,
            'backward': None
        }
        if not isinstance(self.initial_data, (np.ndarray, str)) and self.initial_data is not None:
            state['initial_data'] = {key: self.initial_data[key] for key in ['BR', 'PT', 'TY name', 'TY number', 'Label']}
        return state

    @abstractmethod
    def _set_from_dict(self, state, load_initial_data=True):
        pass

    @contextmanager
    def temporary_chdir(self, new_dir):
        """
            As AUTOCommands does not provide functionality to pass a filepath, this is a workaround to switch paths for loading/saving.
        """
        # Store the current directory
        original_dir = os.getcwd()
        if self._path_name is not None:
            os.chdir(new_dir)
        try:
            yield
        finally:
            # Restore the original directory
            os.chdir(original_dir)
    
    def auto_save(self, auto_suffix):
        if self.continuation:
            for direction in ['forward', 'backward']:
                if self.continuation[direction] is not None:
                    # TODO: I would rather not manually change directory like this, but I can find no option in AUTOcommands to pass a filepath
                    with self.temporary_chdir(self._path_name):
                        ac.save(self.continuation[direction], auto_suffix + '_' + direction)
        self.auto_filename_suffix = auto_suffix

    def auto_load(self, auto_suffix):
        self.continuation = dict()
        for direction in ['forward', 'backward']:
            try:
                # Changing the directory
                # TODO: I would rather not manually change directory like this, but I can find no option in AUTOcommands to pass a filepath
                with self.temporary_chdir(self._path_name):
                    r = ac.loadbd(auto_suffix + '_' + direction)
                self.continuation[direction] = r
            except (FileNotFoundError, OSError, AUTORuntimeError):
                self.continuation[direction] = None

        if self.continuation['forward'] is None and self.continuation['backward'] is None:
            warnings.warn('Files not found. Unable to load data.')
        else:
            self.auto_filename_suffix = auto_suffix

    def save(self, filename=None, auto_filename_suffix=None, **kwargs):
        if auto_filename_suffix is None:
            if self.isfixedpoint:
                self.auto_filename_suffix = "fp_"+str(self.branch_number)
            else:
                self.auto_filename_suffix = "po_" + str(self.branch_number)
            warnings.warn('No AUTO filename suffix set. Using a default one: ' + self.auto_filename_suffix)
        else:
            self.auto_filename_suffix = auto_filename_suffix
        self.auto_save(self.auto_filename_suffix)
        if filename is None:
            if self.isfixedpoint:
                filename = "fp_"+str(self.branch_number)+'.pickle'
            else:
                filename = "po_"+str(self.branch_number)+'.pickle'
            warnings.warn('No pickle filename prefix provided. Using a default one: ' + filename)
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
            if self.continuation['forward'] is not None:
                return self.continuation['forward'].data[0].keys()
            elif self.continuation['backward'] is not None:
                return self.continuation['backward'].data[0].keys()
            else:
                return None
        else:
            return None

    def diagnostics(self):
        if self.continuation:
            s = ''
            if self.continuation['forward'] is not None:
                s += 'Forward\n-----------------------\n'
                for t in self.continuation['forward'].data[0].diagnostics:
                    s += t['Text']
            if self.continuation['backward'] is not None:
                s += 'Backward\n-----------------------\n'
                for t in self.continuation['backward'].data[0].diagnostics:
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
                    if self.continuation['backward'] is not None:
                        s = self.get_solution_by_label(idx)
                        idx = s['PT']
                        ix_map = self._solutions_index_map(direction='backward')
                        if idx is not None:
                            return self.continuation['backward'].data[0].diagnostics[ix_map[idx]]['Text']
                        else:
                            warnings.warn('No point diagnostic to show.')
                            return None
                    else:
                        warnings.warn('No backward branch to show the diagnostic for.')
                        return None
                else:
                    if self.continuation['forward'] is not None:
                        s = self.get_solution_by_label(idx)
                        idx = s['PT']
                        ix_map = self._solutions_index_map(direction='forward')
                        if idx is not None:
                            return self.continuation['forward'].data[0].diagnostics[ix_map[idx]]['Text']
                        else:
                            warnings.warn('No point diagnostic to show.')
                            return None
                    else:
                        warnings.warn('No forward branch to show the diagnostic for.')
                        return None

            if idx >= 0:
                if self.continuation['forward'] is not None:
                    ix_map = self._solutions_index_map(direction='forward')
                    if idx in ix_map:
                        return self.continuation['forward'].data[0].diagnostics[ix_map[idx]]['Text']
                    else:
                        warnings.warn('Point index not found. No point diagnostic to show.')
                        return None
                else:
                    warnings.warn('No forward branch to show the diagnostic for.')
                    return None
            else:
                if self.continuation['backward'] is not None:
                    ix_map = self._solutions_index_map(direction='backward')
                    if -idx in ix_map:
                        return self.continuation['backward'].data[0].diagnostics[ix_map[-idx]]['Text']
                    else:
                        warnings.warn('Point index not found. No point diagnostic to show.')
                        return None
                else:
                    warnings.warn('No backward branch to show the diagnostic for.')
                    return None
        else:
            return None

    def find_solution_index(self, label):
        s = self.get_solution_by_label(label)
        return s['PT']

    def _find_solution_python_auto_index(self, label):
        if self.continuation:
            if label[0] == "-":
                if self.solutions_label['backward'] is not None:
                    label = label[1:]
                    tgt = label[:2]
                    sn = int(label[2:])
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
                    return self._solutions_python_auto_index['backward'][idx]
                else:
                    return None
            else:
                if self.solutions_label['forward'] is not None:
                    tgt = label[:2]
                    sn = int(label[2:])
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
                    return self._solutions_python_auto_index['forward'][idx]
                else:
                    return None
        else:
            return None

    @staticmethod
    def _parse_diagnostic(diag):
        # Extract Point Number from text
        extracted_vals = list()

        lines = diag.strip().splitlines()
        rows = [re.split(r'\s{1,}', line.strip()) for line in lines if line.strip()]
        val_num = 1
        pt_num = -1

        for line in rows:
            if "Multiplier" in line:
                if len(line) == 9:
                    # Assumes a split array of [BR, PT, "Multiplier", multiplier num, FM_re, FM_im, "Abs.", "Val", ab_val]
                    if pt_num == -1:
                        pt_num = int(line[1])
                    else:
                        if pt_num != int(line[1]):
                            raise UserWarning("Cannot find consistent Point Number")

                    if val_num == int(line[3]):
                        # In the case that some numbers are passed as Fortran HUGE(1.0D0), we try two methods to take floats
                        try:
                            re_val = float(line[4])
                        except ValueError:
                            if re.match(r'^-?\d+(\.\d+)?[+-]\d+$', line[4]):
                                re_val = float(re.sub(r'([+-]\d+)$', r'e\1', line[4]))
                            else:
                                raise UserWarning("unexpected form of float")
                        try:
                            im_val = float(line[5])
                        except ValueError:
                            if re.match(r'^-?\d+(\.\d+)?[+-]\d+$', line[5]):
                                im_val = float(re.sub(r'([+-]\d+)$', r'e\1', line[5]))
                            else:
                                raise UserWarning("unexpected form of float")

                        extracted_vals.append([re_val, im_val])
                        val_num += 1

            elif "Eigenvalue" in line:
                if len(line) == 6:
                    # Assumes a split array of [BR, PT, "Multiplier", multiplier num":", lambda_re, lambda_im]
                    if pt_num == -1:
                        pt_num = int(line[1])
                    else:
                        if pt_num != int(line[1]):
                            raise UserWarning("Cannot find consistent Point Number")

                    if val_num == int(line[3][:-1]):
                        # In the case that some numbers are passed as Fortran HUGE(1.0D0), we try two methods to take floats
                        try:
                            re_val = float(line[4])
                        except ValueError:
                            if re.match(r'^-?\d+(\.\d+)?[+-]\d+$', line[4]):
                                re_val = float(re.sub(r'([+-]\d+)$', r'e\1', line[4]))
                            else:
                                raise UserWarning("unexpected form of float")
                        try:
                            im_val = float(line[5])
                        except ValueError:
                            if re.match(r'^-?\d+(\.\d+)?[+-]\d+$', line[5]):
                                im_val = float(re.sub(r'([+-]\d+)$', r'e\1', line[5]))
                            else:
                                raise UserWarning("unexpected form of float")

                        extracted_vals.append([re_val, im_val])
                        val_num += 1
            else:
                continue

        # if pt_num has not been updated, attempt to update it using the first row of data:
        if rows[0][1] == "PT":
            pt_num = int(rows[1][1])

        return extracted_vals, pt_num

    def _solutions_index_map(self, direction='forward'):
        """
            Function creates a map between the `Point number` as found in the solution file, and the `Point number` in the diagnostic d. file.
            It was found that when AUTO cannot converge, it still logs the point and the d. file index then does not correspond with the sol file.
        """
        ix_map = dict()

        if self.continuation[direction] is not None:
            diag_data = self.continuation[direction].data[0].diagnostics.__dict__.copy()['data']

            for i, d in enumerate(diag_data):
                stab_vals, pt_num = self._parse_diagnostic(d['Text'])
                if i < len(diag_data) - 1:
                    dd = diag_data[i + 1]
                    stab_vals_next, pt_num_next = self._parse_diagnostic(dd['Text'])
                    if len(stab_vals) > 0 and pt_num != pt_num_next:
                        ix_map[pt_num] = i
                else:
                    # Catch to compare final entry of diagnostic
                    if len(stab_vals) > 0:
                        ix_map[pt_num] = i
        return ix_map

    @property
    def stability(self):
        if self.continuation:
            s = dict()
            for direction in ['forward', 'backward']:
                if self.continuation[direction] is not None:
                    s[direction] = self.continuation[direction].data[0].stability()
                else:
                    s[direction] = list()
            return s
        else:
            return None

    @property
    def number_of_points(self):
        if self.continuation:
            n = dict()
            for direction in ['forward', 'backward']:
                if self.continuation[direction] is not None:
                    n[direction] = abs(self.continuation[direction].data[0].stability()[-1])
                else:
                    n[direction] = 0
            return n
        else:
            return None

    @property
    def continuation_parameters(self):
        if self.continuation:
            if self.continuation['forward']:
                return self.continuation['forward'].c['ICP']
            elif self.continuation['backward']:
                return self.continuation['backward'].c['ICP']
            else:
                return list()
        else:
            return list()

    @property
    def solutions_index(self):
        d = self.solutions_list_by_direction
        rd = dict()
        if d is not None:
            rd['forward'] = list()
            for sol in d['forward']:
                idx_pt = sol['PT']
                rd['forward'].append(abs(idx_pt))
            rd['backward'] = list()
            for sol in d['backward']:
                idx_pt = sol['PT']
                rd['backward'].append(abs(idx_pt))
            return rd
        else:
            return None


    @property
    def _solutions_python_auto_index(self):
        if self.continuation:
            d = dict()
            if self.continuation['forward'] is not None:
                idx = self.continuation['forward'].data[0].labels.getIndices()
                d['forward'] = idx
            else:
                d['forward'] = list()
            if self.continuation['backward'] is not None:
                idx = self.continuation['backward'].data[0].labels.getIndices()
                d['backward'] = idx
            else:
                d['backward'] = list()
            return d
        else:
            return None

    @property
    def solutions_label(self):
        indices = self._solutions_python_auto_index
        if indices is not None:
            d = dict()
            if indices['forward']:
                idx = indices['forward']
                sl = list()
                for i in idx:
                    sl.append(str(self.continuation['forward'].data[0].labels.by_index[i].keys()).split("'")[1])
                d['forward'] = sl
            else:
                d['forward'] = list()
            if indices['backward']:
                idx = indices['backward']
                sl = list()
                for i in idx:
                    sl.append(str(self.continuation['backward'].data[0].labels.by_index[i].keys()).split("'")[1])
                d['backward'] = sl
            else:
                d['backward'] = list()
            return d
        else:
            return None

    @property
    def number_of_solutions(self):
        sd = dict()
        sols = self.solutions_list_by_direction
        if self.continuation:
            for direction in ['forward', 'backward']:
                sd[direction] = len(sols[direction])
            return sd
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
            if self.continuation['forward'] is not None:
                lab_list = list()
                for lab in self.solutions_label['forward']:
                    if lab not in lab_list:
                        lab_list.append(lab)
                for lab in lab_list:
                    sd['forward'].extend(self.continuation['forward'].getLabel(lab))
            if self.continuation['backward'] is not None:
                lab_list = list()
                for lab in self.solutions_label['backward']:
                    if lab not in lab_list:
                        lab_list.append(lab)
                for lab in lab_list:
                    sd['backward'].extend(self.continuation['backward'].getLabel(lab))
        return sd

    @property
    def solutions_list_by_direction(self):
        indices = self._solutions_python_auto_index
        labels = self.solutions_label
        sd = dict()
        if indices is not None:
            sd['forward'] = list()
            sd['backward'] = list()
            if self.continuation['forward'] is not None:
                for lab, idx in zip(labels['forward'], indices['forward']):
                    if 'solution' in self.continuation['forward'].data[0].labels.by_index[idx][lab]:
                        sol = self.continuation['forward'].data[0].labels.by_index[idx][lab]['solution']
                        sd['forward'].append(sol)
            if self.continuation['backward'] is not None:
                for lab, idx in zip(labels['backward'], indices['backward']):
                    if 'solution' in self.continuation['backward'].data[0].labels.by_index[idx][lab]:
                        sol = self.continuation['backward'].data[0].labels.by_index[idx][lab]['solution']
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

    def get_solution_by_label(self, label):
        if self.continuation:
            if label[0] == '-' and self.continuation['backward'] is not None:
                label = label[1:]
                return self.continuation['backward'].getLabel(label)
            elif label[0] != '-' and self.continuation['forward'] is not None:
                return self.continuation['forward'].getLabel(label)
            else:
                return None
        else:
            return None

    def get_solution_by_index(self, idx):
        if self.continuation:
            sd = self.solutions_list_by_direction
            if idx >= 0:
                for s in sd['forward']:
                    if s['PT'] == idx:
                        return s
                else:
                    warnings.warn(f'Solution for index {idx} not found.')
            elif idx < 0:
                for s in sd['backward']:
                    if s['PT'] == abs(idx):
                        return s
                else:
                    warnings.warn(f'Solution for index {idx} not found.')
            else:
                return None
        else:
            return None

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
            solutions_list = new_solutions_list

        return solutions_list

    def plot_branch_parts(self, variables=(0, 1), ax=None, figsize=(10, 8), markersize=12., plot_kwargs=None, marker_kwargs=None,
                          excluded_labels=('UZ', 'EP'), plot_sol_points=False, variables_name=None, cmap=None):

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

        if plot_sol_points:
            if isinstance(cmap, str):
                cmap = plt.get_cmap(cmap)
            if cmap is None:
                cmap = plt.get_cmap('Reds')

        for direction in ['forward', 'backward']:

            if self.continuation[direction] is not None:

                labels = list()
                for j, coords in enumerate(zip(self.continuation[direction][var1], self.continuation[direction][var2])):
                    lab = self.continuation[direction].data[0].getIndex(j)['TY name']
                    if lab == 'No Label':
                        pass
                    else:
                        labels.append((coords, lab))

                pidx = 0
                for idx in self.stability[direction]:
                    if idx < 0:
                        ls = '-'
                    else:
                        ls = '--'
                    plot_kwargs['linestyle'] = ls
                    lines_list = ax.plot(self.continuation[direction][var1][pidx:abs(idx)], self.continuation[direction][var2][pidx:abs(idx)], **plot_kwargs)
                    if plot_sol_points:
                        # generate points and colours for plotting
                        x_vals = self.solutions_parameters(parameters=var1)[0]
                        p_min, p_max = np.min(x_vals), np.max(x_vals)
                        y_vals = np.empty_like(x_vals)
                        point_color = list()
                        for i, x in enumerate(x_vals):
                            ix = np.argmin(np.abs(self.continuation[direction][var1][pidx:abs(idx)] - x))
                            y_vals[i] = (self.continuation[direction][var2][pidx:abs(idx)][ix])
                            point_color.append(cmap((x - p_min) / (p_max - p_min)))

                        ax.scatter(x_vals, y_vals, c=point_color)

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

    def plot_branch_parts_3D(self, variables=(0, 1, 2), ax=None, figsize=(10, 8), markersize=12., plot_kwargs=None, marker_kwargs=None,
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

        for direction in ['forward', 'backward']:

            if self.continuation[direction] is not None:

                labels = list()
                for j, coords in enumerate(zip(self.continuation[direction][var1], self.continuation[direction][var2], self.continuation[direction][var3])):
                    lab = self.continuation[direction].data[0].getIndex(j)['TY name']
                    if lab == 'No Label':
                        pass
                    else:
                        labels.append((coords, lab))

                pidx = 0
                for idx in self.stability[direction]:
                    if idx < 0:
                        ls = '-'
                    else:
                        ls = '--'
                    plot_kwargs['linestyle'] = ls
                    lines_list = ax.plot(self.continuation[direction][var1][pidx:abs(idx)], self.continuation[direction][var2][pidx:abs(idx)],
                                         self.continuation[direction][var3][pidx:abs(idx)], **plot_kwargs)
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
                       linewidth=None, color_solutions=False, plot_kwargs=None, labels=None, indices=None, parameter=None, value=None,
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

        # Colouring solutions dependant on parameter value
        if 'cmap' in plot_kwargs:
            cmap = plot_kwargs['cmap']
            if isinstance(cmap, str):
                cmap = plt.get_cmap(cmap)
            plot_kwargs.pop('cmap')
        else:
            cmap = plt.get_cmap('Blues')
        if color_solutions and (parameter is not None):
            p_vals = self.solutions_parameters(parameters=parameter)
            p_min, p_max = np.min(p_vals), np.max(p_vals)

            if value is None:
                value = p_vals

        solutions_list = self.get_filtered_solutions_list(labels, indices, parameter, value, None, tol)

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
            if color_solutions and (parameter is not None):
                plot_kwargs['color'] = cmap((sol.PAR[parameter] - p_min) / (p_max - p_min))
            line_list = ax.plot(x, y, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth,
                                **plot_kwargs)

            if not color_solutions:
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

        solutions_list = self.get_filtered_solutions_list(labels, indices, parameter, value, None, tol)

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
        if self.continuation['forward'] is not None:
            summary_str += "Forward\n"
            summary_str += "=======\n"
            summary_str += self.continuation['forward'].summary() + "\n\n"

        if self.continuation['backward'] is not None:
            summary_str += "Backward\n"
            summary_str += "========\n"
            summary_str += self.continuation['backward'].summary()

        return summary_str

    def print_summary(self):

        if self.continuation:
            s = self.summary()
            print(s)

    def check_for_repetitions(self, parameters, tol=2.e-2, return_parameters=False, return_non_repeating_solutions=False,
                              return_repeating_solutions=False, forward=None):

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        repeating_solution = dict()
        repeating_solution['backward'] = list()
        repeating_solution['forward'] = list()
        if forward:
            repeating_solution['backward'] = None
        elif not forward:
            repeating_solution['forward'] = None

        forward_dict = {'forward': True, 'backward': False}

        sol = dict()
        idx_dict = dict()
        ridx_dict = dict()
        for direction in ['forward', 'backward']:
            if repeating_solution[direction] is not None:
                idx_dict[direction] = list()
                ridx_dict[direction] = list()
                ridx_list = list()
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
                            ridx_list.append(i)
                            break

                for i in range(ssol.shape[1]):
                    if i not in ridx_list:
                        repeating_solution[direction].append(False)
                        idx_dict[direction].append(i)

                    else:
                        repeating_solution[direction].append(True)
                        ridx_dict[direction].append(i)

        if forward is None:
            if repeating_solution['backward'] is not None:
                res = [repeating_solution['backward']]
            else:
                res = list([])
            if repeating_solution['forward'] is not None:
                res[0] += repeating_solution['forward']
        elif forward:
            if repeating_solution['forward'] is not None:
                res = [repeating_solution['forward']]
            else:
                res = list([])
        else:
            if repeating_solution['backward'] is not None:
                res = [repeating_solution['backward']]
            else:
                res = list([])

        if return_parameters:
            if forward is None:
                params = np.concatenate((sol['backward'][:, ridx_dict['backward']], sol['foward'][:, ridx_dict['forward']]), axis=-1)
            elif forward:
                params = sol['forward'][:, idx_dict['forward']]
            else:
                params = sol['backward'][:, idx_dict['backward']]

            res.append(params)

        for i, t in enumerate([return_non_repeating_solutions, return_repeating_solutions]):
            if t:
                tarr = dict()
                for direction in ['backward', 'forward']:
                    if repeating_solution[direction] is not None:
                        tarr[direction] = np.array(repeating_solution[direction])
                        if i == 0:
                            tarr[direction] = ~tarr[direction]

                all_sols = self.solutions_list_by_direction
                sols = list()
                if forward is None:
                    for direction in ['backward', 'forward']:
                        if repeating_solution[direction] is not None:
                            dir_sols = np.empty(len(all_sols[direction]), dtype=object)
                            dir_sols[:] = all_sols[direction]
                            sel_sols = dir_sols[tarr[direction]]
                            sols.extend(sel_sols)
                else:
                    if forward:
                        direction = 'forward'
                    else:
                        direction = 'backward'
                    if repeating_solution[direction] is not None:
                        dir_sols = np.empty(len(all_sols[direction]), dtype=object)
                        dir_sols[:] = all_sols[direction]
                        sel_sols = dir_sols[tarr[direction]]
                        sols.extend(sel_sols)
                res.append(sols)

        if len(res) == 1:
            return res[0]
        else:
            return res


def _sort_arrays(sol, npar, tol):

    if npar == 1:
        srt = np.squeeze(np.argsort(sol))
    else:
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

