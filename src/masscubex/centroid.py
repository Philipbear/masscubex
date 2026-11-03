"""Profile MS1 mass detection with local background and physical m/z moments."""

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def centroid_profile(mz, intensity, baseline, threshold, min_points=3, gap_factor=4.0):
    """Return centroid m/z, raw apex intensity and apex noise threshold.

    Resolve local maxima (including flat tops), split at intervening valleys,
    and stop at noise crossings or discontinuous sampling. The centroid is the
    first moment of the piecewise-linear, background-subtracted mass peak.
    At least ``min_points`` contiguous samples above threshold are required;
    isolated one-bin electronic spikes and truncated edge maxima are excluded.
    This general-purpose detector is not a vendor peak-picking replacement.
    """
    n = mz.size
    peaks = np.empty(n, dtype=np.int64)
    starts = np.empty(n, dtype=np.int64)
    ends = np.empty(n, dtype=np.int64)
    count = 0
    i = 1
    while i < n - 1:
        if intensity[i] <= intensity[i - 1]:
            i += 1
            continue
        right = i
        while right + 1 < n and intensity[right + 1] == intensity[i]:
            right += 1
        if right < n - 1 and intensity[right] > intensity[right + 1]:
            apex = (i + right) // 2
            peaks[count] = apex
            starts[count], ends[count] = i, right
            count += 1
        i = right + 1
    valleys = np.empty(count + 1, dtype=np.int64)
    valleys[0] = 0
    valleys[count] = n - 1
    for k in range(1, count):
        a, b = ends[k - 1], starts[k]
        valley = a
        for j in range(a + 1, b + 1):
            if intensity[j] < intensity[valley]:
                valley = j
        valleys[k] = valley
    out_mz = np.empty(count)
    out_intensity = np.empty(count)
    out_threshold = np.empty(count)
    kept = 0
    for k in range(count):
        apex = peaks[k]
        if intensity[apex] <= threshold[apex]:
            continue
        spacing = min(mz[starts[k]] - mz[starts[k] - 1], mz[ends[k] + 1] - mz[ends[k]])
        max_gap = gap_factor * spacing
        left = right = apex
        while left > valleys[k]:
            if intensity[left - 1] <= threshold[left - 1] or mz[left] - mz[left - 1] > max_gap:
                break
            left -= 1
        while right < valleys[k + 1]:
            if intensity[right + 1] <= threshold[right + 1] or mz[right + 1] - mz[right] > max_gap:
                break
            right += 1
        if right - left + 1 < min_points:
            continue
        # Integrate in coordinates relative to the apex to preserve m/z accuracy.
        origin = mz[apex]
        area = moment = 0.0
        for j in range(left, right):
            x0, x1 = mz[j] - origin, mz[j + 1] - origin
            y0 = max(0.0, intensity[j] - baseline[j])
            y1 = max(0.0, intensity[j + 1] - baseline[j + 1])
            dx = x1 - x0
            area += dx * (y0 + y1) / 2
            moment += dx * (x0 * (2 * y0 + y1) + x1 * (y0 + 2 * y1)) / 6
        if area <= 0:
            continue
        out_mz[kept] = origin + moment / area
        out_intensity[kept] = intensity[apex]
        out_threshold[kept] = threshold[apex]
        kept += 1
    return out_mz[:kept], out_intensity[:kept], out_threshold[:kept]
