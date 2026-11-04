"""Independent Gaussian diagnostic for calibration-kernel tests."""

import numpy as np


def reference_gaussian_similarity(times, values):
    if len(values) < 5 or values.max() <= 0:
        return 0.0
    apex = int(np.argmax(values))
    y = values / values[apex]
    before = np.flatnonzero(y[:apex] < 0.5)
    after = np.flatnonzero(y[apex:] < 0.5)
    left = int(before[-1]) if len(before) else 0
    right = apex + int(after[0]) if len(after) else len(y) - 1
    sigma = max((times[right] - times[left]) / 2.355, 1e-6)
    model = np.exp(-0.5 * ((times - times[apex]) / sigma) ** 2)
    if np.std(y) == 0 or np.std(model) == 0:
        return 0.0
    return float(np.clip(np.corrcoef(y, model)[0, 1], 0, 1))
