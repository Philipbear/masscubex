"""Explicit, immutable settings shared safely by independent file workers."""

from dataclasses import asdict, dataclass, fields
import math

from .adducts import Adduct, DEFAULT_ADDUCT_IDS, selected_adducts


@dataclass(frozen=True)
class Config:
    scan_cache_mb: float = 64.0
    auto_tune: bool = True
    tuning_min_peaks: int = 20
    tuning_max_peaks: int = 100
    tuning_search_ppm: float = 100.0
    tuning_search_da: float = 0.05
    tuning_peak_limit: int = 128
    tuning_block_scans: int = 512
    tuning_block_seconds: float = 120.0
    tuning_locked: tuple[str, ...] = ()
    mz_tolerance_ppm: float = 10.0
    mz_tolerance_da: float = 0.002
    noise_mode: str = "adaptive"
    noise_estimator: str = "low_tail"
    noise_layout: str = "bins"
    noise_fraction: float = 0.10
    noise_window_da: float = 100.0
    noise_sigma: float = 3.0
    noise_clip_sigma: float = 2.0
    noise_iterations: int = 2
    noise_min_neighbors: int = 5
    min_intensity: float = 0.0
    min_scans: int = 5
    max_gap_scans: int = 2
    max_gap_seconds: float = 15.0
    segment: bool = True
    smoothing_sigma: float = 0.6
    smoothing_seconds: float = 0.0
    peak_prominence: float = 0.1
    peak_distance: int = 3
    peak_distance_seconds: float = 0.0
    assume_centroid: bool = False
    profile_mode: str = "centroid"
    profile_min_points: int = 3
    profile_gap_factor: float = 4.0
    quant_method: str = "peak_height"
    associate_ms2: bool = True
    ms2_mz_tolerance_ppm: float = 20.0
    ms2_mz_tolerance_da: float = 0.01
    ms2_rt_tolerance: float = 2.0
    ms2_max_isolation_width: float = 0.0
    ms2_seeded_features: bool = True
    ms2_gapfill_rt_window: float = 5.0
    ms2_seed_rt_tolerance: float = 5.0
    gapfill_mz_tolerance_ppm: float = 10.0
    gapfill_mz_tolerance_da: float = 0.002
    gapfill_noise_fraction: float = 0.1
    precursor_max_age_seconds: float = 15.0
    gapfill_max_gap_seconds: float = 15.0
    best_ms2: str = "fragment_intensity"
    group_ions: bool = True
    group_isf: bool = True
    adduct_ids: tuple[int, ...] = DEFAULT_ADDUCT_IDS
    custom_adducts: tuple[Adduct, ...] = ()
    group_mz_tolerance_ppm: float = 10.0
    group_mz_tolerance_da: float = 0.002
    group_rt_tolerance: float = 3.0
    group_min_scans: int = 5
    group_min_correlation: float = 0.9
    group_min_overlap: float = 0.5
    isf_min_relative_intensity: float = 0.05
    isf_max_peaks: int = 100

    def __post_init__(self):
        integer_minima = {
            "tuning_min_peaks": 20,
            "tuning_max_peaks": 20,
            "tuning_peak_limit": 32,
            "tuning_block_scans": 16,
            "noise_iterations": 0,
            "noise_min_neighbors": 2,
            "min_scans": 2,
            "max_gap_scans": 0,
            "peak_distance": 1,
            "profile_min_points": 3,
            "group_min_scans": 3,
            "isf_max_peaks": 1,
        }
        for name, minimum in integer_minima.items():
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        booleans = {
            "auto_tune",
            "segment",
            "assume_centroid",
            "associate_ms2",
            "ms2_seeded_features",
            "group_ions",
            "group_isf",
        }
        for name in booleans:
            if type(getattr(self, name)) is not bool:
                raise ValueError(f"{name} must be a boolean")
        if self.noise_mode not in {"adaptive", "none"}:
            raise ValueError("noise_mode must be 'adaptive' or 'none'")
        if self.noise_estimator not in {"low_tail", "clipped"}:
            raise ValueError("noise_estimator must be 'low_tail' or 'clipped'")
        if self.noise_layout not in {"bins", "sliding"}:
            raise ValueError("noise_layout must be 'bins' or 'sliding'")
        if self.profile_mode not in {"centroid", "reject"}:
            raise ValueError("profile_mode must be 'centroid' or 'reject'")
        if self.quant_method not in {"peak_height", "peak_area"}:
            raise ValueError("quant_method must be peak_height or peak_area")
        if self.best_ms2 not in {"fragment_intensity", "precursor_intensity", "apex"}:
            raise ValueError("best_ms2 must be fragment_intensity, precursor_intensity or apex")
        if not isinstance(self.custom_adducts, tuple):
            raise ValueError("custom_adducts must be a tuple")
        selected_adducts(self.adduct_ids, self.custom_adducts)
        if self.tuning_max_peaks < self.tuning_min_peaks:
            raise ValueError("tuning_max_peaks must be >= tuning_min_peaks")
        if not isinstance(self.tuning_locked, tuple) or any(
            key not in TUNABLE_FIELDS for key in self.tuning_locked
        ):
            raise ValueError("tuning_locked must contain valid tunable parameter names")
        nonnegative = {
            "scan_cache_mb",
            "min_intensity",
            "mz_tolerance_da",
            "mz_tolerance_ppm",
            "ms2_mz_tolerance_da",
            "ms2_mz_tolerance_ppm",
            "ms2_rt_tolerance",
            "ms2_max_isolation_width",
            "group_mz_tolerance_da",
            "group_mz_tolerance_ppm",
            "group_rt_tolerance",
            "isf_min_relative_intensity",
            "smoothing_seconds",
            "peak_distance_seconds",
            "gapfill_mz_tolerance_da",
            "gapfill_mz_tolerance_ppm",
            "gapfill_noise_fraction",
        }
        excluded = (
            set(integer_minima)
            | booleans
            | {
                "noise_mode",
                "noise_estimator",
                "noise_layout",
                "profile_mode",
                "quant_method",
                "best_ms2",
                "adduct_ids",
                "custom_adducts",
                "tuning_locked",
            }
        )
        for field in fields(self):
            if field.name in excluded:
                continue
            value = getattr(self, field.name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{field.name} must be a finite number")
            if not math.isfinite(value) or (value < 0 if field.name in nonnegative else value <= 0):
                raise ValueError(f"Invalid {field.name}: {value}")
        if self.noise_sigma != 3.0:
            raise ValueError("noise_sigma is fixed at 3; adjust noise_fraction to tune filtering")
        if self.mz_tolerance_da == 0 and self.mz_tolerance_ppm == 0:
            raise ValueError("At least one m/z tolerance must be positive")
        for prefix in ("ms2", "group"):
            if (
                getattr(self, prefix + "_mz_tolerance_da") == 0
                and getattr(self, prefix + "_mz_tolerance_ppm") == 0
            ):
                raise ValueError(f"At least one {prefix} m/z tolerance must be positive")
        if self.gapfill_mz_tolerance_da == 0 and self.gapfill_mz_tolerance_ppm == 0:
            raise ValueError("At least one gapfill m/z tolerance must be positive")
        for name in (
            "noise_fraction",
            "gapfill_noise_fraction",
            "group_min_correlation",
            "group_min_overlap",
            "isf_min_relative_intensity",
        ):
            if getattr(self, name) > 1:
                raise ValueError(f"{name} must be <= 1")
        if self.peak_prominence > 1:
            raise ValueError("peak_prominence must be <= 1")

    def to_dict(self):
        result = asdict(self)
        result["adduct_ids"] = list(self.adduct_ids)
        result["custom_adducts"] = [asdict(ion) for ion in self.custom_adducts]
        result["tuning_locked"] = list(self.tuning_locked)
        return result

    @classmethod
    def from_dict(cls, values):
        if not isinstance(values, dict):
            raise ValueError("Configuration must be a JSON object")
        unknown = set(values) - {field.name for field in fields(cls)}
        if unknown:
            raise ValueError(f"Unknown configuration keys: {', '.join(sorted(unknown))}")
        values = dict(values)
        try:
            if "adduct_ids" in values:
                values["adduct_ids"] = tuple(values["adduct_ids"])
            if "tuning_locked" in values:
                values["tuning_locked"] = tuple(values["tuning_locked"])
            if "custom_adducts" in values:
                values["custom_adducts"] = tuple(Adduct(**ion) for ion in values["custom_adducts"])
            return cls(**values)
        except (TypeError, KeyError) as error:
            raise ValueError(f"Invalid configuration: {error}") from error


TUNABLE_FIELDS = frozenset(
    {
        "mz_tolerance_ppm",
        "mz_tolerance_da",
        "ms2_mz_tolerance_ppm",
        "ms2_mz_tolerance_da",
        "gapfill_mz_tolerance_ppm",
        "gapfill_mz_tolerance_da",
        "max_gap_seconds",
        "smoothing_seconds",
        "peak_distance_seconds",
        "ms2_rt_tolerance",
        "ms2_gapfill_rt_window",
        "ms2_seed_rt_tolerance",
        "precursor_max_age_seconds",
        "gapfill_max_gap_seconds",
        "group_rt_tolerance",
    }
)
