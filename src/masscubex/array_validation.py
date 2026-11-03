"""One-pass checks for the common, already normalized spectrum-array case."""

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def normalized_arrays(mz, intensity, strict_mass_order, positive_intensity):
    """Call after checking dimensions/lengths; never mutate either input array."""
    previous = 0.0
    for index in range(len(mz)):
        mass, value = mz[index], intensity[index]
        if not np.isfinite(mass) or not np.isfinite(value) or mass <= 0:
            return False
        if value < 0 or (positive_intensity and value == 0):
            return False
        if mass < previous or (strict_mass_order and mass == previous):
            return False
        previous = mass
    return True
