# Data-derived PD drive topography, variants (a) GPe->STN potentiation only and (b) beta depth only (the other hook
# uniform: synapse 2.0 / depth M). The uniform controls are identical across variants (ctl_*, run by
# drive_topography_c.py) and are not repeated. See sweeps/drive_topography_common.py.
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drive_topography_common import *  # noqa: F401,F403,E402
from drive_topography_common import CONTROLS, D, pd_condition  # noqa: E402

conditions = [pd_condition(n, v) for v in ('a', 'b') for n in D.CONDITIONS if n not in CONTROLS]
