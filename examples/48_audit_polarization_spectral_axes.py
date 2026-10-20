"""Audit GHOST polarization products and spectral-channel structure with CASA."""

from pathlib import Path
import json

import numpy as np
from casatools import table

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "data/ghost_f1_5min/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
OUT = ROOT / "data/ghost_matrix/extension_smoke_tests/ms_axis_audit.json"

tb = table()
tb.open(str(MS / "POLARIZATION"))
try:
    corr_type = np.asarray(tb.getcol("CORR_TYPE"), int).ravel().tolist()
    num_corr = np.asarray(tb.getcol("NUM_CORR"), int).ravel().tolist()
finally:
    tb.close()
tb.open(str(MS / "SPECTRAL_WINDOW"))
try:
    num_chan = np.asarray(tb.getcol("NUM_CHAN"), int).ravel().tolist()
    chan_freq = np.asarray(tb.getcol("CHAN_FREQ"), float)
    chan_width = np.asarray(tb.getcol("CHAN_WIDTH"), float)
finally:
    tb.close()

names = {1:"I", 2:"Q", 3:"U", 4:"V", 5:"RR", 6:"RL", 7:"LR", 8:"LL",
         9:"XX", 10:"XY", 11:"YX", 12:"YY"}
result = {"measurement_set":str(MS.relative_to(ROOT)), "corr_type_codes":corr_type,
          "corr_type_names":[names.get(v, str(v)) for v in corr_type],
          "num_corr":num_corr, "num_chan":num_chan,
          "chan_frequency_hz_shape":list(chan_freq.shape),
          "chan_width_hz_shape":list(chan_width.shape),
          "has_cross_hands":any(v in corr_type for v in (6, 7, 10, 11))}
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result, indent=2))
