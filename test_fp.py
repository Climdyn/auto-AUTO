
# temporary test file - to be removed later

from auto2.parsers.config import ConfigParser
from auto2.continuation.fixed_points import FixedPointContinuation
from auto2.continuation.periodic_orbits import PeriodicOrbitContinuation
import numpy as np

c = ConfigParser('c.qgs_land-atmosphere_auto')

fp = FixedPointContinuation('qgs_land-atmosphere_auto', c)

ic = np.zeros(c.ndim)

fp.make_continuation(ic, ICP=['C_go1'], PAR={'C_go1': 0., 2: 0., 3: 0.085, 4: 0.02})

# fp.auto_save('fp1')

s = fp.get_filtered_solutions_list(labels='BP')[1]

fp2 = FixedPointContinuation('qgs_land-atmosphere_auto', c)
fp2.make_continuation(s, ISW=-1)

s = fp.get_filtered_solutions_list(labels='HB')[0]

hp = PeriodicOrbitContinuation('qgs_land-atmosphere_auto', c)

hp.make_continuation(s, max_bp=None, ICP=['C_go1', 'T'], IPS=2)

s2 = fp.get_filtered_solutions_list(labels='HB')[1]

hp2 = PeriodicOrbitContinuation('qgs_land-atmosphere_auto', c)

hp2.make_continuation(s2, max_bp=None, ICP=['C_go1', 'T'], IPS=2, NMX=385)
