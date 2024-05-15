
# temporary test file - to be removed later

from auto2.parsers.config import ConfigParser
from auto2.continuation.fixed_points import FixedPointContinuation
import numpy as np

c = ConfigParser('c.qgs_land-atmosphere_auto')

fp = FixedPointContinuation('qgs_land-atmosphere_auto', c)

ic = np.zeros(c.ndim)

fp.make_continuation(ic, ICP=['C_go1'], PAR={1: 0., 2: 0., 3: 0.085, 4: 0.02})

# fp.auto_save('fp1')
