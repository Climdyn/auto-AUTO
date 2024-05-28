
# temporary test file - to be removed later

import numpy as np
from numba import njit
from scipy.optimize import root

from auto2.diagrams.bifurcations import BifurcationDiagram

from qgs.params.params import QgParams
from qgs.functions.tendencies import create_tendencies

model_parameters = QgParams({'phi0_npi': np.deg2rad(50.)/np.pi, 'n':1.3 }, dynamic_T=False)

model_parameters.set_atmospheric_channel_fourier_modes(2, 2)
model_parameters.set_ground_channel_fourier_modes(2, 2)

# Changing (increasing) the orography depth
model_parameters.ground_params.set_orography(0.2, 1)
# Setting the parameters of the heat transfer from the soil
model_parameters.gotemperature_params.set_params({'gamma': 1.6e7, 'T0': 300})
model_parameters.atemperature_params.set_params({ 'hlambda':10, 'T0': 290})
# Setting atmospheric parameters
model_parameters.atmospheric_params.set_params({'sigma': 0.2, 'kd': 0.085, 'kdp': 0.02})

# Setting insolation 
model_parameters.gotemperature_params.set_params({})

C_g = 300
model_parameters.atemperature_params.set_insolation(0.4*C_g , 0)

model_parameters.gotemperature_params.set_insolation(C_g , 0)

f, Df = create_tendencies(model_parameters)

@njit
def fnt(x):
    return f(0., x)


nsearch = 1000

# Start on random initial conditions
ic = 2 * (np.random.rand(nsearch, model_parameters.ndim) - 0.5) * 0.05

eps = 1.e-6
fixed_points = dict()

sol_idx = 1
for i in range(nsearch):
    sol = root(fnt, ic[i, :])
    if sol.success:
        for idx in fixed_points:
            if np.linalg.norm(fixed_points[idx] - sol.x) < eps:
                break
        else:
            fixed_points[sol_idx] = sol.x
            sol_idx+=1

par = {'C_go1': 300.}

initial_points = list()

for p in fixed_points:
    initial_points.append({'parameters': par, 'initial_data': fixed_points[p]})

print('Found ' + str(len(fixed_points)) + ' fixed points. Computing bifurcation diagram.')

b = BifurcationDiagram('qgs_land-atmosphere_auto')

b.compute_fixed_points_diagram(initial_points, extra_comparison_parameters=['psi_a_1', 'psi_a_2', 'psi_a_5'], comparison_tol=[1.e-3, 1.e-4, 1.e-4, 1.e-4],
                               ICP=['C_go1'], PAR={'C_go1': 300., 2: 300. * 0.4, 3: 0.085, 4: 0.02})
# b.compute_fixed_points_diagram(initial_points, extra_comparison_parameters=['psi_a_2'], comparison_tol=[2.e-2, 2.e-3],
#                                ICP=['C_go1'], PAR={'C_go1': 300., 2: 300. * 0.4, 3: 0.085, 4: 0.02})
