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
        sys.path.append(auto_directory + '/python/auto')
        sys.path.append(auto_directory + '/python')
except KeyError:
    warnings.warn('Unable to find auto directory environment variable.')

import AUTOCommands as ac
import runAUTO as ra

# TODO: Add solutions plot
# TODO: Add diagnostics for given point
# TODO: Add stability info for given point
# TODO: Allow starting from solution
# TODO: Define as subclass of base class


class FixedPointContinuation(object):

    def __init__(self, model_name, config_object):

        self.config_object = config_object
        self.model_name = model_name
        self.continuation = list()
        self.branch_number = None

    def make_continuation(self, initial_data, store_name="", only_forward=False, **continuation_kwargs):
        runner = ra.runAUTO()
        ac.load(self.model_name, runner=runner)
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

    def plot_branches(self, variables=(0, 1), ax=None, figsize=(10, 8), markersize=12., plot_kwargs=None, marker_kwargs=None,
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

        keys = self.available_variables

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

        for branches in [0, 1]:

            if self.continuation[branches] is not None:

                labels = list()
                for j, coords in enumerate(zip(self.continuation[branches][var1], self.continuation[branches][var2])):
                    lab = self.continuation[branches].data[0].getIndex(j)['TY name']
                    if lab == 'No Label':
                        pass
                    else:
                        labels.append((coords, lab))

                pidx = 0
                for idx in self.stability[branches]:
                    if idx < 0:
                        ls = '-'
                    else:
                        ls = '--'
                    plot_kwargs['ls'] = ls
                    lines_list = ax.plot(self.continuation[branches][var1][pidx:abs(idx)], self.continuation[branches][var2][pidx:abs(idx)], **plot_kwargs)
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

