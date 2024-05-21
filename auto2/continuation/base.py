from abc import ABC, abstractmethod
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


class Continuation(ABC):

    def __init__(self, model_name, config_object):
        self.config_object = config_object
        self.model_name = model_name
        self.continuation = list()
        self.branch_number = None
        self.initial_data = None

        # plots default behaviours
        self._default_marker = None
        self._default_markersize = None
        self._default_linestyle = None
        self._default_linewidth = None

    @abstractmethod
    def make_continuation(self, initial_data, store_name="", only_forward=False, **continuation_kwargs):
        pass

    def auto_save(self, store_name):
        if self.continuation:
            ac.save(self.continuation[0], store_name + '_forward')
            if self.continuation[1] is not None:
                ac.save(self.continuation[1], store_name + '_backward')

    def auto_load(self, store_name):
        self.continuation = list()
        r = ac.loadbd(store_name + '_forward')
        self.continuation.append(r)
        try:
            r = ac.loadbd(store_name + '_backward')
            self.continuation.append(r)
        except:
            self.continuation.append(None)


    @property
    def available_variables(self):
        if self.continuation:
            return self.continuation[0].data[0].keys()
        else:
            return None

    def diagnostics(self):
        if self.continuation:
            s = 'Forward\n-----------------------\n'
            for t in self.continuation[0].data[0].diagnostics:
                s += t['Text']
            if self.continuation[1] is not None:
                s = 'Backward\n-----------------------\n'
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
                    warnings.warn('No backward branch to show the diagnostic for.')
                    return None

        if idx >= 0:
            return self.continuation[0].data[0].diagnostics[idx]
        else:
            if self.continuation[1] is not None:
                return self.continuation[1].data[0].diagnostics[-idx]['Text']
            else:
                warnings.warn('No backward branch to show the diagnostic for.')
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

    @property
    def stability(self):
        if self.continuation:
            s = list()
            s.append(self.continuation[0].data[0].stability())
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
            n.append(self.continuation[0].data[0].stability()[1])
            if self.continuation[1] is not None:
                n.append(self.continuation[1].data[0].stability()[1])
            else:
                n.append(0)
            return n
        else:
            return None

    @property
    def solutions_index(self):
        if self.continuation:
            d = dict()
            idx = self.continuation[0].data[0].labels.getIndices()
            d['forward'] = idx
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
            idx = indices['forward']
            sl = list()
            for i in idx:
                sl.append(str(self.continuation[0].data[0].labels.by_index[i].keys()).split("'")[1])
            d['forward'] = sl
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
            sl.append(self.continuation[0].data[0].getLabels()[-1])
            if self.continuation[1] is not None:
                sl.append(self.continuation[1].data[0].getLabels()[-1])
            else:
                sl.append(0)
            return sl
        else:
            return None

    @property
    def full_solutions_list(self):
        sd = self.solutions_list_by_direction
        if sd['backward']:
            sl = sd['backward'][1:]
        else:
            sl = list()
        sl.extend(sd['forward'])
        return sl

    @property
    def solutions_list_by_direction(self):
        sd = dict()
        if self.solutions_label is not None:
            sd['forward'] = list()
            sd['backward'] = list()
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

    def solutions_parameters(self, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD'), forward=None):
        if not isinstance(parameters, (tuple, list)):
            parameters = [parameters]
        sl = self.get_filtered_solutions_list(labels=solutions_types, forward=forward)
        params = dict()
        for param in parameters:
            par_list = list()
            for s in sl:
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
            for sol in solutions_list:
                for val in values.T:
                    for i, param in enumerate(parameters):
                        if abs(sol[param] - val[i]) > tol[i]:
                            break
                    else:
                        new_solutions_list.append(sol)
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

    def same_solutions_as(self, other, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD'), tol=2.e-2):

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
            dif = ssol - osol
            return np.all(np.abs(dif).T < tol)

    def solutions_in(self, other, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD'), tol=2.e-2, return_parameters=False, return_solutions=False):
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

    def solutions_part_of(self, other, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD'), tol=2.e-2, return_parameters=False, return_solutions=False, forward=None):

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        ssol = self.solutions_parameters(parameters, solutions_types, forward=forward)
        osol = other.solutions_parameters(parameters, solutions_types)

        npar = ssol.shape[0]
        if isinstance(tol, float):
            tol = [tol] * npar

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        ssol = _sort_arrays(ssol, npar, tol)
        osol = _sort_arrays(osol, npar, tol)

        idx_list = list()
        for i in range(ssol.shape[1]):
            for j in range(osol.shape[1]):
                dif = ssol[:, i] - osol[:, j]
                if np.all(np.abs(dif) < tol):
                    idx_list.append(i)

        res = [len(idx_list) > 1]
        if return_parameters:
            res.append(ssol[:, idx_list])
        if return_solutions:
            res.append(self.get_filtered_solutions_list(parameters=parameters, values=ssol[:, idx_list], tol=tol))

        if len(res) == 1:
            return res[0]
        else:
            return res

    def branch_possibly_cross(self, other, parameters, solutions_types=('HB', 'BP', 'UZ', 'PD'), tol=2.e-2, return_parameters=False, return_solutions=False, forward=None):

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        ssol = self.solutions_parameters(parameters, solutions_types, forward=forward)
        osol = other.solutions_parameters(parameters, solutions_types)

        npar = ssol.shape[0]
        if isinstance(tol, float):
            tol = [tol] * npar

        if isinstance(tol, (list, tuple)):
            tol = np.array(tol)

        ssol = _sort_arrays(ssol, npar, tol)
        osol = _sort_arrays(osol, npar, tol)

        idx_list = list()
        for i in range(ssol.shape[1]):
            for j in range(osol.shape[1]):
                dif = ssol[:, i] - osol[:, j]
                if np.all(np.abs(dif) < tol):
                    idx_list.append(i)

        res = [len(idx_list) == 1]
        if return_parameters:
            res.append(ssol[:, idx_list])
        if return_solutions:
            res.append(self.get_filtered_solutions_list(parameters=parameters, values=ssol[:, idx_list], tol=tol)[0])

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
