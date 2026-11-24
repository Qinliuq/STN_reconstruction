# Data-derived PD drive topography, variant (c): both hooks scaled. Uniform controls (shared by all variants) and the
# new healthy references. See sweeps/drive_topography_common.py.
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drive_topography_common import *  # noqa: F401,F403,E402
from drive_topography_common import CONTROLS, D, healthy_conditions, pd_condition  # noqa: E402

conditions = [pd_condition(n, 'c') for n in D.CONDITIONS] + healthy_conditions()
