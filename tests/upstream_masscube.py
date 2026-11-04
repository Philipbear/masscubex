"""Unmodified function bodies from huaxuyu/masscube, Apache-2.0.

Copyright 2024-2026 Huaxu Yu. Commit 356748ede691862b9a87f67812db5e728a703b9d.
Test-only source snapshot: adapters supply spectra and the selected ion catalog.
Do not edit/format these function bodies when changing MassCubeX.
"""

# ruff: noqa
# fmt: off
from __future__ import annotations
from copy import deepcopy
import bisect
import numpy as np
from numpy import dot
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

def tqdm(values):
    return values

POS_ADDUCTS = {}
NEG_ADDUCTS = {}


# Source: src/masscube/feature_evaluation.py:14
def calculate_gaussian_similarity(x, y, len_tol=5):
    """
    calculate the gaussian similarity using Pearson correlation coefficient.

    Parameters
    ----------
    x: numpy array
        Retention times.
    y: numpy array
        MS1 signal intensities.
    len_tol: int
        Minimum length of the peak to calculate the similarity score. 
        If the length of the peak is less than this value, return 0.
    
    Returns
    -------
    float
        similarity score
    """

    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if not isinstance(y, np.ndarray):
        y = np.array(y)

    # if the length of the peak is less than the tolerance, return 0
    if len(x) < len_tol or np.max(y) <= 0:
        return 0.0

    # Estimate the parameters of the Gaussian function
    # a: amplitude, b: mean, c: standard deviation
    a = np.max(y)
    idx = np.argmax(y)
    b = x[idx]

    # Compute approximate FWHM using half-height
    right_idx = np.where(y[idx:] < a / 2)[0]
    left_idx = np.where(y[:idx] < a / 2)[0]

    c1 = (x[right_idx[0] + idx] - b) if len(right_idx) > 0 else (x[-1] - b)
    c2 = (b - x[left_idx[-1]]) if len(left_idx) > 0 else (b - x[0])

    fwhm = c1 + c2
    c = fwhm / 2.355 if fwhm > 0 else 1e-6  # Avoid division by zero

    # Generate fitted Gaussian curve
    y_fit = a * np.exp(-0.5 * ((x - b) / c) ** 2)

    # Fast correlation using dot product
    y = y.astype(float)
    y_fit = y_fit.astype(float)
    y -= y.mean()
    y_fit -= y_fit.mean()
    similarity = dot(y, y_fit) / (np.linalg.norm(y) * np.linalg.norm(y_fit) + 1e-12)
    
    return max(0, similarity)

# Source: src/masscube/feature_evaluation.py:72
def calculate_noise_score(
        y: np.ndarray, 
        rel_int_tol: float = 0.05, 
        len_tol: int = 5
    ) -> float:
    """
    Noise score ranges from 0 to infinity and reflects the fluctuation of the peak. 
    A higher noise score indicates a noisier peak.

    Parameters
    ----------
    y: np.ndarray
        Intensity values of the peak.
    rel_int_tol: float
        Relative intensity threshold (filter low signals).
    len_tol: int
        Minimum length of the peak to calculate the noise score. 
        If the length of the peak is less than this value, return 0.
    
    Returns
    -------
    float
        Noise score (0 = no noise, higher = noiser)
    """

    y = np.asarray(y, dtype=float).copy()

    if len(y) < len_tol or np.max(y) <= 0:
        return 0.0

    y[y < rel_int_tol * np.max(y)] = 0  # Filter low signals

    total_variation = np.abs(np.diff(y)).sum()
    total_variation += y[0] + y[-1]

    return total_variation / (2 * y.max()) - 1

# Source: src/masscube/feature_evaluation.py:110
def calculate_asymmetry_factor(y: np.ndarray, threshold_ratio: float = 0.1) -> float:
    """
    Calculate peak asymmetry at a given threshold level (default: 10%).

    Parameters
    ----------
    y : np.ndarray
        Intensity values
    threshold_ratio : float
        Fraction of peak height to define tails (default = 0.1 = 10%)

    Returns
    -------
    float
        Asymmetry factor:
        - ~1 = symmetric
        - >1 = right-skewed
        - <1 = left-skewed
        - np.inf = undefined
    """

    y = np.array(y, dtype=float) if not isinstance(y, np.ndarray) else y

    if len(y) < 5 or np.max(y) <= 0:
        return np.inf

    idx = np.argmax(y)

    # Define threshold
    threshold = y[idx] * threshold_ratio

    # Find last point before peak below threshold (left tail)
    left_candidates = np.where(y[:idx] < threshold)[0]
    left_idx = left_candidates[-1] if len(left_candidates) else 0

    # Find first point after peak below threshold (right tail)
    right_candidates = np.where(y[idx:] < threshold)[0]
    right_idx = (right_candidates[0] + idx) if len(right_candidates) else (len(y) - 1)

    # Avoid division by zero
    if idx == left_idx:
        return 99.0

    return (right_idx - idx) / (idx - left_idx)

# Source: src/masscube/feature_detection.py:20
class Feature:
    """
    A class to store a feature characterized by a unique pair of m/z and retention time.
    """

    def __init__(self):

        # chromatographic peak
        self.rt_seq = []                     # retention time sequence
        self.signals = []                    # signal sequence (m/z, intensity)
        self.scan_idx_seq = []               # scan index sequence
        self.ms2_seq = []                    # MS2 spectra
        self.gap_counter = 0                 # count the number of consecutive zeros in the end of the peak
        self.peak_shape = None               # peak shape ([[rt, intensity], ...]) 
        self.peak_edges = None               # retention time edges for the peak baseline: (left_rt, right_rt)

        # summary
        self.id = None                       # feature id
        self.feature_group_id = None         # peak group id
        self.mz = None                       # m/z
        self.rt = None                       # retention time
        self.scan_idx = None                 # scan index of the peak apex
        self.peak_height = None              # peak height
        self.peak_area = None                # peak area
        self.top_average = None              # average of the highest three intensities
        self.ms2 = None                      # representative MS2 spectrum (default: highest total intensity)
        self.length = 0                      # number of valid scans in the feature
        self.gaussian_similarity = 0.0       # Gaussian similarity
        self.noise_score = 0.0               # noise score
        self.asymmetry_factor = 0.0          # asymmetry factor
        self.sse = 0.0                       # squared error to the smoothed curve
        self.is_segmented = False            # whether the feature is segmented from a larger feature
        self.is_isotope = False              # whether the feature is an isotope
        self.charge_state = 1                # charge state of the feature
        self.isotope_signals = None          # isotope signals [[mz, intensity], ...]
        self.is_in_source_fragment = False   # whether the feature is an in-source fragment
        self.adduct_type = None              # adduct type
        self.isotope_state = "M+0"           # isotope state, e.g., M+0, M+1, M+2, etc.
        self.scan_scan_cor = None            # scan-scan correlation to the most abundant ion in the same feature group

        self.annotation_algorithm = None     # annotation algorithm. Not used now.
        self.search_mode = None              # 'identity search', 'fuzzy search', or 'mzrt_search'
        self.similarity = None               # similarity score (0-1)
        self.annotation = None               # name of annotated compound
        self.formula = None                  # molecular formula
        self.matched_peak_number = None      # number of matched peaks
        self.smiles = None                   # SMILES
        self.inchikey = None                 # InChIKey
        self.matched_precursor_mz = None     # matched precursor m/z
        self.matched_ms2 = None              # matched ms2 spectra
        self.matched_adduct_type = None      # matched adduct type


    def extend(self, rt, signal, scan_idx):
        """
        Function to extend the chromatographic peak by providing the retention time, 
        signal and scan index.

        Parameters
        ----------
        rt: float
            Retention time.
        signal: list or numpy array
            [mz, intensity] of the signal.
        scan_idx: int
            Scan index.
        """

        self.rt_seq.append(rt)
        self.signals.append(signal)
        self.scan_idx_seq.append(scan_idx)


    def get_mz_error(self):
        """
        Function to calculate the 3*sigma error of the feature's m/z.

        Returns
        -------
        float
            The 3*sigma error of the feature's m/z.
        """

        return 3 * np.std(self.signals[:, 0])
    

    def get_rt_error(self):
        """
        Function to calculate the 3*sigma error of the feature's retention time.

        Returns
        -------
        float
            The 3*sigma error of the feature's retention time.
        """

        return 3 * np.std(self.rt_seq)


    def finalize(self, ph=True, pa=True, ta=True, g_score=True, n_score=True, a_score=True):
        """
        Finalize the feature by calculating the summary statistics, including peak height, 
        peak area, top average, Gaussian similarity, noise score, and asymmetry factor.

        Parameters
        ----------
        ph: bool
            Whether to calculate the peak height.
        pa: bool
            Whether to calculate the peak area.
        ta: bool
            Whether to calculate the top average.
        g_score: bool
            Whether to calculate the Gaussian similarity.
        n_score: bool
            Whether to calculate the noise score.
        a_score: bool
            Whether to calculate the asymmetry factor.
        """

        self.signals = np.array(self.signals, dtype=np.float32)
        original_rt_seq = self.rt_seq
        first,last = _trim_signals(self.signals)
        if last > first:
            left_edge_idx = first - 1 if first > 0 else first
            right_edge_idx = last if last < len(original_rt_seq) else last - 1
            self.peak_edges = (original_rt_seq[left_edge_idx], original_rt_seq[right_edge_idx])
        self.signals = self.signals[first:last]
        self.rt_seq = self.rt_seq[first:last]
        self.scan_idx_seq = self.scan_idx_seq[first:last]
        
        apx = np.argmax(self.signals[:, 1])
        self.mz = self.signals[apx, 0]
        self.rt = self.rt_seq[apx]
        self.scan_idx = self.scan_idx_seq[apx]
        self.length = np.sum(self.signals[:, 1] > 0)
        self.peak_shape = np.array([[self.rt_seq[i], self.signals[i, 1]] for i in range(len(self.signals))])

        if ph:
            self.peak_height = int(self.signals[apx, 1])
        if pa:
            self.peak_area = int(np.trapz(y=self.signals[:, 1], x=self.rt_seq) * 60)
        if ta:
            self.top_average = int(np.mean(np.sort(self.signals[:, 1])[-3:]))
        if n_score:
            self.noise_score = calculate_noise_score(self.signals[:, 1])
        if g_score:
            self.gaussian_similarity = calculate_gaussian_similarity(self.rt_seq, self.signals[:, 1])
        if a_score:
            self.asymmetry_factor = calculate_asymmetry_factor(self.signals[:, 1])


    def subset(self, start, end, finalize=True):
        """
        Keep the subset of the feature by providing the start and end positions. Note that the 
        summary statistics will be recalculated in this function by default.

        Parameters
        ----------
        start: int
            The start position.
        end: int
            The end position. The data point at the end position is not included.
        finalize: bool
            Whether to finalize the feature after subsetting. If True, the summary statistics will be recalculated.
        """
        
        self.rt_seq = self.rt_seq[start:end]
        self.signals = self.signals[start:end]
        self.scan_idx_seq = self.scan_idx_seq[start:end]
        self.ms2_seq = [ms2 for ms2 in self.ms2_seq if ms2.raw_file_id > self.scan_idx_seq[0] and ms2.raw_file_id < self.scan_idx_seq[-1]]

        if finalize:
            self.finalize(pa=False, ta=False, g_score=False, a_score=False)

# Source: src/masscube/feature_detection.py:280
def segment_feature(feature, method="gf-prominence", length_tol=5, noise_tol=3):
    """
    Function to segment a feature into multiple features based on the edge detection.

    Parameters
    ----------
    sigma: float
        The sigma value for Gaussian filter.
    prominence_ratio: float
        The prominence ratio for finding peaks. prom = np.max(y)*prominence_ratio
    distance: int
        The minimum distance between peaks.
    peak_height_tol: float
        The peak height tolerance for segmentation.
    length_tol: int
        The length tolerance for segmentation.
    sse_tol: float
        The squared error tolerance for segmentation.

    Returns
    -------
    segmented_features: list
        A list of segmented features.
    """

    # if the number of valid data points is less than length_tol, skip segmentation
    arr = feature.signals[:,1]
    dp = np.sum(arr > 0)
    if dp < length_tol:
        return [feature]
    
    noise_score = calculate_noise_score(arr)
    if noise_score > noise_tol:
        return [feature]

    if method == "gf-prominence":
        peaks = gaussian_filter_and_prominence_method(arr)

    if len(peaks) < 2:
        return [feature]
    
    positions = [0]
    for i in range(len(peaks)-1):
        # get the index of the lowest intensity point between peaks[i] and peaks[i+1]
        lowest_int_idx = peaks[i] + np.argmin(arr[peaks[i]:peaks[i+1]])
        positions.append(lowest_int_idx)
    positions.append(len(feature.signals)-1)
    
    segmented_features = []
    for i in range(len(positions)-1):
        tmp = deepcopy(feature)
        tmp.subset(start=positions[i], end=positions[i+1]+1)
        tmp.is_segmented = True
        segmented_features.append(tmp)

    return segmented_features

# Source: src/masscube/feature_detection.py:338
def gaussian_filter_and_prominence_method(arr, sigma=0.6, prominence_ratio=0.03, distance=3):
    """
    Peak segmentation using Gaussian filter coupled with peak prominence.

    Parameters
    ----------
    arr: numpy array
        1D array of intensity values.
    sigma: float
        The sigma value for Gaussian filter.
    prominence_ratio: float
        The prominence ratio for finding peaks. prom = np.max(y)*prominence_ratio
    distance: int
        The minimum distance between peaks.

    Returns
    -------
    peaks: numpy array
        Indices of detected peaks.
    """
    
    # add zero to the front and the end of the signal to facilitate the edge detection
    a = np.zeros(len(arr)+2)
    a[1:-1] = arr
    ss = gaussian_filter1d(a, sigma=sigma)
    prominence = np.max(ss)*prominence_ratio

    peaks, _ = find_peaks(ss, prominence=prominence, distance=distance)

    peaks = peaks - 1   # correct the index
    
    # the resulting peaks should have peak height larger than the baseline
    peaks = peaks[arr[peaks] > np.median(arr)]

    return peaks

# Source: src/masscube/feature_detection.py:423
def _trim_signals(signals: np.ndarray) -> tuple:
    """
    Return a (first, last) half-open slice that removes leading/trailing zeros
    in the intensity column (assumed to be column index 1).

    Parameters
    ----------
    signals : np.ndarray
        A 2D numpy array of signals shaped like [[mz, intensity], ...].

    Returns
    -------
    first : int
        Index of the first non-zero intensity.
    last : int
        One past the index of the last non-zero intensity (half-open slice).
        If all intensities are zero, returns (0, 0).
    """

    non_zero_indices = np.where(signals[:, 1] != 0)[0]
    if non_zero_indices.size == 0:
        return 0, 0
    else:
        return non_zero_indices[0], non_zero_indices[-1] + 1

# Source: src/masscube/feature_grouping.py:239
def group_features_single_file(d: MSData, rt_tol: float = None) -> None:
    """
    Group isotopes, adducts, and in-source fragments from a single file 
    based on the m/z, retention time, MS/MS, scan-to-scan correlation and relative intensity.
    The feature_group_id attribute of the Feature objects in the MSData object will be updated.

    Parameters
    ----------
    d: MSData object
        Contains the raw data and the features to be grouped.
    rt_tol: float
        The retention time tolerance for grouping features (0.05 min by default).
    """

    if rt_tol is None:
        rt_tol = d.params.rt_tol_feature_grouping

    # check if features have been detected
    if d.features is None or len(d.features) == 0:
        raise ValueError("No features found in the MSData object for feature grouping.")

    # sort features by peak height from high to low
    d.features.sort(key=lambda x: x.peak_height, reverse=True)
    
    # initialization
    mz_arr = np.array([f.mz for f in d.features])
    rt_arr = np.array([f.rt for f in d.features])
    is_grouped = np.zeros(len(d.features), dtype=bool)
    feature_group_id = 1
    ion_mode = d.metadata.ion_mode
    if ion_mode is None:
        raise ValueError("ion_mode metadata is required for single-file feature grouping.")
    default_adduct = "[M+H]+" if ion_mode.lower() == "positive" else "[M-H]-"

    # find isotopes, adducts and in-source fragments for each feature
    for i, f in enumerate(tqdm(d.features)):

        # skip if already grouped
        if is_grouped[i]:
            continue
        
        # mark as grouped
        is_grouped[i] = True
        f.feature_group_id = feature_group_id
        f.isotope_state = "M+0"
        # check if adduct type is assigned by MS/MS library search
        # if not, assign default adduct type based on ionization mode
        if f.adduct_type is None:
            f.adduct_type = default_adduct

        peak_main = f.signals[:,1]
        rt_range  = [f.rt_seq[0] - 1e-6, f.rt_seq[-1] + 1e-6]

        # extract isotope signals from the apex scan
        f.isotope_signals = find_isotope_signals(mz=f.mz, signals=d.scans[f.scan_idx].signals, 
                                                 mz_tol=d.params.mz_tol_feature_grouping, 
                                                 rel_int_limit=d.params.isotope_rel_int_limit)
        
        # masks
        # intensity mask is not needed since features have been sorted by peak height
        rt_mask = np.abs(rt_arr - f.rt) < rt_tol
        base_mask = rt_mask & (~is_grouped)

        # generate all possible ions for grouping
        ion_dict = enumerate_source_ions(mz=f.mz, ms2=f.ms2.signals if f.ms2 is not None else None, 
                                         adduct=f.adduct_type, ion_mode=ion_mode)

        for key, mz_val in ion_dict.items():
            mz_mask = np.abs(mz_arr - mz_val) < d.params.mz_tol_feature_grouping
            mask = base_mask & mz_mask
            matched_num = np.sum(mask)

            # one isotope mz can be matched to multiple features because of high-resolution MS
            if key.startswith("M+"):
                is_grouped[mask] = True
                for j in np.where(mask)[0]:
                    d.features[j].feature_group_id = f.feature_group_id
                    d.features[j].is_isotope = True
                    d.features[j].adduct_type = f.adduct_type
                    d.features[j].isotope_state = key
                    d.features[j].isotope_signals = f.isotope_signals
                continue
            
            # only one matched is allowed for adducts and in-source fragments
            if matched_num > 0:
                # prioritize the one with the lowest RT difference
                v = np.where(mask)[0]
                vi = v[np.argmin(np.abs(rt_arr[v] - f.rt))]
                _, eic, _ = d.get_eic_data(mz_arr[vi], f.rt, rt_range=rt_range)
                scan_scan_cor = scan_to_scan_cor_intensity(peak_main, eic[:, 1])
                
                if scan_scan_cor > d.params.scan_scan_cor_tol:
                    is_grouped[vi] = True
                    d.features[vi].feature_group_id = f.feature_group_id
                    d.features[vi].scan_scan_cor = scan_scan_cor
                    # for adducts:
                    if not key.startswith("ISF"):
                        d.features[vi].adduct_type = key
                    # for in-source fragments:
                    else:
                        d.features[vi].is_in_source_fragment = True
                        d.features[vi].adduct_type = f.adduct_type
                    # find isotopes for this ion
                    d.features[vi].isotope_signals = find_isotope_signals(d.features[vi].mz,
                                                                          signals=d.scans[d.features[vi].scan_idx].signals,
                                                                          mz_tol=d.params.mz_tol_feature_grouping,
                                                                          rel_int_limit=d.params.isotope_rel_int_limit)
                    # retrieve the isotope signals from the feature list
                    iso_deep_dic = enumerate_source_ions(d.features[vi].mz, ms2=None, adduct=None, ion_mode=None, isotope_only=True)
                    for key_deep, mz_val_deep in iso_deep_dic.items():
                        mz_mask_deep = np.abs(mz_arr - mz_val_deep) < d.params.mz_tol_feature_grouping
                        # use the latest grouping status
                        mask_deep = rt_mask & (~is_grouped) & mz_mask_deep
                        for j in np.where(mask_deep)[0]:
                            # need to refine by relative intensity
                            if d.features[j].peak_height < d.features[vi].peak_height * d.params.isotope_rel_int_limit:
                                is_grouped[j] = True
                                d.features[j].feature_group_id = feature_group_id
                                d.features[j].is_isotope = True
                                d.features[j].adduct_type = d.features[vi].adduct_type
                                d.features[j].isotope_state = key_deep
                                d.features[j].isotope_signals = d.features[vi].isotope_signals

        feature_group_id += 1

# Source: src/masscube/feature_grouping.py:365
def enumerate_source_ions(mz: float, ms2: np.array, adduct: str, ion_mode: str, iso_num: int = 4, isotope_only: bool = False) -> dict:
    """
    Calculate all possible isotopes, adducts, and in-source fragments m/z values for a given m/z.
    The calculation does not compute the isotopes of the adducts and in-source fragments.

    Parameters
    ----------
    mz: float
        The m/z value of the feature.
    ms2: numpy.array
        The MS2 signals as [[m/z, intensity], ...]
    adduct: str
        The adduct form of the feature.
    ion_mode: str
        The ionization mode. "positive" or "negative".
    iso_num: int
        The number of isotopes to consider after the monoisotopic peak.
    isotope_only: bool
        Whether to only calculate isotopes. Default is False.

    Returns
    -------
    dict
        A dictionary with keys as adduct/in-source fragment forms and values as the corresponding m/z values.
    """

    ion_dict = {}
    
    # isotopes
    for i in range(1, iso_num):
        ion_dict[f'M+{i}'] = mz + 1 + (i-1) * 1.003355
    
    if isotope_only:
        return ion_dict

    # adducts
    if ion_mode.lower() == "positive":
        adduct_dic = POS_ADDUCTS
    elif ion_mode.lower() == "negative":
        adduct_dic = NEG_ADDUCTS
    
    if adduct in adduct_dic.keys():
        a = adduct_dic[adduct]
        base_mz = (mz - a.mass_shift) * a.charge / a.mol_multiplier
        for key, val in adduct_dic.items():     
            if key != adduct and val.considered:
                ion_dict[key] = (base_mz * val.mol_multiplier + val.mass_shift) / val.charge

    # in-source fragments
    if ms2 is not None:
        if type(ms2) == str:
            ms2 = extract_signals_from_string(ms2)
        for i, p in enumerate(ms2):
            ion_dict[f'ISF_{i+1}'] = p[0]
    
    return ion_dict

# Source: src/masscube/feature_grouping.py:423
def find_isotope_signals(mz: float, signals: np.array, mz_tol: float = 0.01, 
                         num: int = 5, rel_int_limit: float = 0.7) -> np.array:
    """
    Find isotopes from MS1 signals. Only single charged isotopes are considered.

    Parameters
    ----------
    mz: float
        The m/z value of the feature.
    signals: numpy.array
        The MS1 signals as [[m/z, intensity], ...]
    mz_tol: float
        The m/z tolerance to find isotopes.
    num: int
        The maximum number of isotopes to be found.
    rel_int_limit: float
        Isotope's intensity cannot be higher than rel_int_limit * intensity of the monoisotopic peak.

    Returns
    -------
    numpy.array
        The m/z and intensity of the isotopes.
    """

    mzs = signals[:, 0]
    intens = signals[:, 1]

    idx = np.abs(mzs - mz).argmin()
    limit = intens[idx] * rel_int_limit

    targets = mz + np.arange(num, dtype=mzs.dtype)
    within = np.any(np.abs(mzs[None, :] - targets[:, None]) < mz_tol, axis=0)
    good_int = intens < limit

    keep = within & good_int
    
    return signals[keep]

# Source: src/masscube/feature_grouping.py:467
def scan_to_scan_cor_intensity(a: np.array, b: np.array) -> float:
    """
    Calculate the scan-to-scan correlation (Pearson correlation) between two intensity arrays.

    Parameters
    ----------
    a: np.array
        Intensity array of the first m/z
    b: np.array
        Intensity array of the second m/z
    
    Returns
    -------
    float
        The scan-to-scan correlation (Pearson correlation) between the two intensity arrays.
    """

    v = (a>0) & (b>0)

    # if the commonly detected points are less than 3, ignore the correlation calculation and return 1.0
    if np.sum(v) < 3:
        return 1.0

    return np.corrcoef(a[v], b[v])[0, 1]

# Source: src/masscube/feature_detection.py
def detect_features(d):
    """
    Detect features in the MS data.

    Parameters
    ----------
    d: MSData object
        An MSData object that contains the MS data.

    Returns
    -------
    final_features: list
        A list of detected features.
    """

    # A list to store the rois in progress
    features = []
    # A list for the finally detected ROIs
    final_features = []

    # Initiate a set of rois using the first MS1 scan
    s = d.scans[d.ms1_idx_arr[0]]    # The first scan

    for i in range(len(s.signals)):
        feature = Feature()
        feature.extend(rt=s.time, signal=s.signals[i], scan_idx=d.ms1_idx_arr[0])
        features.append(feature)

    # Loop over all MS1 scans
    for ms1_idx in d.ms1_idx_arr[1:]:
        s = d.scans[ms1_idx]                                  # The current MS1 scan
        avlb_signals = np.ones(len(s.signals), dtype=bool)    # available signals to assign to features
        to_be_moved = []                                      # features to be moved to final_features
        
        for i, feature in enumerate(features):
            min_idx = None
            if len(s.signals) > 0:
                min_idx = _find_closest_index_ordered(
                    array=s.signals[:, 0],
                    target=feature.signals[-1][0],
                    tol=d.params.mz_tol_ms1,
                )
            if min_idx is not None and avlb_signals[min_idx]:
                feature.extend(rt=s.time, signal=s.signals[min_idx], scan_idx=ms1_idx)
                feature.gap_counter = 0
                avlb_signals[min_idx] = False
            else:
                feature.extend(rt=s.time, signal=[feature.signals[-1][0], 0], scan_idx=ms1_idx)
                feature.gap_counter = feature.gap_counter + 1
                if feature.gap_counter > d.params.feature_gap_tol:
                    to_be_moved.append(i)

        # Move the features that have not been visited for a long time to final_features
        for i in to_be_moved[::-1]:
            final_features.append(features.pop(i))
        
        # Create new rois for the remaining signals
        for i, signal in enumerate(s.signals):
            if avlb_signals[i]:
                feature = Feature()
                feature.extend(rt=s.time, signal=signal, scan_idx=ms1_idx)
                features.append(feature)

        features.sort(key=lambda x: x.signals[-1][1], reverse=True)

    # Move all features to final_features
    for feature in features:
        final_features.append(feature)
    
    # finalize features
    for feature in final_features:
        feature.finalize(pa=False, g_score=False, a_score=False)
    
    # sort by m/z
    final_features.sort(key=lambda x: x.mz)

    return final_features

# Source: src/masscube/feature_detection.py
def _find_closest_index_ordered(array, target, tol=0.01):
    """
    Function to find the index of the closest value in an ordered array.

    Parameters
    ----------
    array: list or numpy array
        An ordered array.
    target: float
        The target value.
    tol: float
        The tolerance for the closest value.

    Returns
    -------
    idx: int
        The index of the closest value.
    """

    idx = bisect.bisect_left(array, target)
    
    if idx == 0:
        if array[idx] - target < tol:
            return 0
        else:
            return None
    if idx == len(array):
        if target - array[idx - 1] < tol:
            return len(array) - 1
        else:
            return None
    
    before = array[idx - 1]
    after = array[idx]
    
    if after - target < target - before and after - target < tol:
        return idx
    elif after - target > target - before and target - before < tol:
        return idx - 1
    else:
        return None

# Source: src/masscube/utils_functions.py
def centroid_signals(signals, mz_tol=0.002):
    """
    Function to centroid signals in a mass spectrum.

    Parameters
    ----------
    signals: numpy array
        MS signals for a scan as 2D numpy array in float32, organized as [[m/z, intensity], ...].
    mz_tol: float
        m/z tolerance for centroiding. Default is 0.002 Da.

    Returns
    -------
    signals: numpy array
        Centroided signals.
    """

    if mz_tol is None or len(signals) <= 1:
        return signals

    # Assume signals sorted by m/z
    signals = signals[signals[:, 0].argsort()]
    mz = signals[:, 0]
    intensity = signals[:, 1]

    diff = np.diff(mz)
    group_starts = np.where(diff >= mz_tol)[0] + 1
    group_boundaries = np.r_[0, group_starts, len(signals)]

    sum_intensity = np.add.reduceat(intensity, group_boundaries[:-1])
    weighted_mz = np.add.reduceat(mz * intensity, group_boundaries[:-1]) / sum_intensity

    merged_signals = np.column_stack((weighted_mz, sum_intensity)).astype(np.float32)
    
    return merged_signals

# Source: src/masscube/raw_data_utils/core.py
def find_best_ms2(ms2_list):
    """
    Function to find the best MS2 spectrum for a list of MS2 spectra.
    """

    if len(ms2_list) > 0:
        total_ints = [np.sum(ms2.signals[:,1]) for ms2 in ms2_list]
        if np.max(total_ints) == 0:
            return ms2_list[0]
        else:
            return ms2_list[max(range(len(total_ints)), key=total_ints.__getitem__)]
    else:
        return None

# Source: src/masscube/raw_data_utils/core.py
def _set_feature_peak_edges(self, feature):
        """
        Set feature baseline edges to the adjacent MS1 scan times when available.
        """

        if len(feature.rt_seq) == 0 or len(feature.scan_idx_seq) == 0 or len(self.ms1_idx_arr) == 0:
            return

        left_edge = feature.rt_seq[0]
        right_edge = feature.rt_seq[-1]

        left_pos = np.searchsorted(self.ms1_idx_arr, feature.scan_idx_seq[0])
        if left_pos < len(self.ms1_idx_arr) and self.ms1_idx_arr[left_pos] == feature.scan_idx_seq[0]:
            if left_pos > 0:
                left_edge = self.scans[self.ms1_idx_arr[left_pos - 1]].time

        right_pos = np.searchsorted(self.ms1_idx_arr, feature.scan_idx_seq[-1])
        if right_pos < len(self.ms1_idx_arr) and self.ms1_idx_arr[right_pos] == feature.scan_idx_seq[-1]:
            if right_pos + 1 < len(self.ms1_idx_arr):
                right_edge = self.scans[self.ms1_idx_arr[right_pos + 1]].time

        feature.peak_edges = (left_edge, right_edge)

# Source: src/masscube/raw_data_utils/core.py
def allocate_ms2_to_features(self, mz_tol: float = 0.1) -> None:
        """
        Assign MS2 scans to features by precursor m/z and retention time.

        The deliberately broad default m/z window tolerates raw files with
        inaccurate precursor-ion m/z metadata. Candidate features within that
        window are disambiguated by peak height.
        
        If multiple features are matched, the MS2 scan is assigned to the
        feature with the highest peak height. Each MS2 scan is assigned to at
        most one feature.

        Parameters
        ----------
        mz_tol : float, default=0.1
            m/z tolerance for matching MS2 precursor m/z to feature m/z.
        """

        if len(self.features) == 0 or len(self.ms2_idx_arr) == 0:
            return

        feature_mz_arr = np.asarray(self.feature_mz_arr, dtype=np.float64)

        sorted_idx = np.argsort(feature_mz_arr)
        sorted_feature_mz_arr = feature_mz_arr[sorted_idx]

        feature_rt_start_arr = np.array(
            [
                feature.peak_edges[0]
                if feature.peak_edges is not None
                else feature.rt_seq[0]
                for feature in self.features
            ],
            dtype=np.float32,
        )
        feature_rt_end_arr = np.array(
            [
                feature.peak_edges[1]
                if feature.peak_edges is not None
                else feature.rt_seq[-1]
                for feature in self.features
            ],
            dtype=np.float32,
        )
        feature_height_arr = np.array([feature.peak_height for feature in self.features], dtype=np.float32)

        for scan_idx in self.ms2_idx_arr:
            ms2 = self.scans[scan_idx]

            left = np.searchsorted(sorted_feature_mz_arr, ms2.precursor_mz - mz_tol, side="left")
            right = np.searchsorted(sorted_feature_mz_arr, ms2.precursor_mz + mz_tol, side="right")
            
            if left == right:
                continue

            candidate_idx = sorted_idx[left:right]
            
            mz_mask = np.abs(feature_mz_arr[candidate_idx] - ms2.precursor_mz) < mz_tol
            rt_mask = (
                (feature_rt_start_arr[candidate_idx] <= ms2.time)
                & (ms2.time <= feature_rt_end_arr[candidate_idx])
            )
            matched_idx = candidate_idx[mz_mask & rt_mask]
            if len(matched_idx) == 0:
                continue

            best_idx = matched_idx[np.argmax(feature_height_arr[matched_idx])]
            ms2.ms2_allocated = True
            self.features[best_idx].ms2_seq.append(ms2)

# Source: src/masscube/raw_data_utils/core.py
def get_eic_data(self, target_mz, target_rt=None, mz_tol=0.005, rt_tol=0.3, rt_range=None):
        if rt_range is None:
            if target_rt is None:
                rt0, rt1 = 0.0, np.inf
            else:
                rt0, rt1 = target_rt - rt_tol, target_rt + rt_tol
        else:
            rt0, rt1 = rt_range[0], rt_range[1]

        times = self.ms1_time_arr          # (n_ms1,) float32/float64, sorted
        ms1_idx = self.ms1_idx_arr         # (n_ms1,) int32, aligned to times
        scans = self.scans                 # local bind

        # RT window -> contiguous slice via binary search (much faster than mask+where)
        left = np.searchsorted(times, rt0, side="left")
        right = np.searchsorted(times, rt1, side="right")
        if right <= left:
            return _EMPTY_F32, _EMPTY_SIG, _EMPTY_I32

        eic_time_arr = times[left:right]
        eic_scan_idx_arr = ms1_idx[left:right]
        n = eic_scan_idx_arr.size

        # allocate outputs as two 1D arrays (faster than (n,2) then column ops)
        eic_mz = np.full(n, np.nan, dtype=np.float32)
        eic_int = np.zeros(n, dtype=np.float32)

        mz0 = float(target_mz)
        lo = mz0 - mz_tol
        hi = mz0 + mz_tol

        for out_i, scan_i in enumerate(eic_scan_idx_arr):
            sig = scans[int(scan_i)].signals
            if sig is None or sig.shape[0] == 0:
                continue

            mzs = sig[:, 0]  # sorted ascending
            # find m/z window indices in O(log n_peaks)
            l = np.searchsorted(mzs, lo, side="left")
            r = np.searchsorted(mzs, hi, side="right")
            if r <= l:
                continue

            ints = sig[:, 1]
            j = l + int(np.argmax(ints[l:r]))
            eic_mz[out_i] = mzs[j]
            eic_int[out_i] = ints[j]

        eic_signals = np.column_stack((eic_mz, eic_int)).astype(np.float32, copy=False)
        return eic_time_arr, eic_signals, eic_scan_idx_arr
