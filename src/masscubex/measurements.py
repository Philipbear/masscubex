"""Specialize repeated measurements without changing their NumPy arithmetic.

The general scientific reference remains in upstream.py. These paths reuse a
known apex and the fixed, unweighted two-vector form of NumPy's covariance and
correlation calculations. Reductions and BLAS operations retain their order;
there is no fastmath, approximate Pearson formula or changed peak-shape model.
"""

# Covariance/correlation arithmetic adapted from NumPy 1.26 (BSD-3-Clause).
# Copyright (c) 2005-2023, NumPy Developers. See licenses/NUMPY-BSD-3-Clause.txt.

import numpy as np


def integrated_area(minutes, values):
    """Upstream integer area: the 1-D np.trapz expression, then minutes to seconds."""
    # Keep the float32 intensity addition before the float64 time product.
    # NumPy's sum order and the final Python int conversion are intentional.
    differences = np.diff(np.asarray(minutes))
    total = (differences * (values[1:] + values[:-1]) / 2.0).sum(axis=-1)
    return float(int(total * 60))


def noise_score_at_apex(values, apex):
    """MassCube's trace-noise score, reusing the unchanged maximum after its 5% filter."""
    if len(values) < 5:
        return 0.0
    maximum = np.float64(values[apex])
    if maximum <= 0:
        return 0.0
    # The reference promotes to float64 and then copies that temporary again.
    # One owned float64 array suffices; no caller data are modified.
    y = values.astype(float)
    y[y < 0.05 * maximum] = 0
    variation = np.abs(np.diff(y)).sum()
    variation += y[0] + y[-1]
    return variation / (2 * maximum) - 1


def gaussian_similarity(times, values, apex):
    """MassCube's Gaussian similarity with the already known first maximum."""
    if len(times) < 5:
        return 0.0
    amplitude = values[apex]
    if amplitude <= 0:
        return 0.0
    x = np.asarray(times)
    center = x[apex]
    right = np.where(values[apex:] < amplitude / 2)[0]
    left = np.where(values[:apex] < amplitude / 2)[0]
    right_width = x[right[0] + apex] - center if len(right) else x[-1] - center
    left_width = center - x[left[-1]] if len(left) else center - x[0]
    fwhm = right_width + left_width
    width = fwhm / 2.355 if fwhm > 0 else 1e-6
    fitted = amplitude * np.exp(-0.5 * ((x - center) / width) ** 2)
    observed = values.astype(float)
    fitted = fitted.astype(float, copy=False)
    observed -= observed.mean()
    fitted -= fitted.mean()
    similarity = np.dot(observed, fitted) / (
        np.linalg.norm(observed) * np.linalg.norm(fitted) + 1e-12
    )
    return max(0, similarity)


def pearson_correlation(x, y):
    """np.corrcoef(x, y)[0, 1] for two real vectors with at least three values.

    Keep NumPy's row means, matrix dot product, sample-covariance scaling,
    sequential standard-deviation divisions and clipping exactly. Only the
    general API's argument validation and weight handling are omitted.
    """
    if np.iscomplexobj(x) or np.iscomplexobj(y):
        return np.corrcoef(x, y)[0, 1]
    samples = np.empty((2, len(x)), dtype=np.result_type(x, y, np.float64))
    samples[0], samples[1] = x, y
    samples -= samples.mean(axis=1)[:, None]
    covariance = np.dot(samples, samples.T.conj())
    covariance *= np.true_divide(1, len(x) - 1)
    deviation = np.sqrt(np.diag(covariance).real)
    covariance /= deviation[:, None]
    covariance /= deviation[None, :]
    np.clip(covariance.real, -1, 1, out=covariance.real)
    return covariance[0, 1]
