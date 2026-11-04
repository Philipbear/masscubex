"""Single-pass offline mzML/mzXML readers with optional MS2 metadata/peaks.

No byte-offset index, whole-file scan list, or downloaded CV is needed. The
reader follows PSI mzML 1.1 and mzXML 2/3 binary conventions. Unsupported
encodings fail explicitly instead of producing an apparently valid table.
"""

from dataclasses import dataclass
from pathlib import Path
import base64
import gzip
import math
import re
import zlib

from lxml import etree
import numpy as np

from .array_validation import normalized_arrays
from .formats import file_format


@dataclass(slots=True)
class Scan:
    rt: float  # seconds
    mz: np.ndarray
    intensity: np.ndarray
    polarity: str = "unknown"
    native_id: str = ""
    centroided: bool = True
    spectrum_index: int = -1
    acquisition_range: tuple[float, float] | None = None
    invalid_acquisition_range: bool = False


@dataclass(slots=True)
class MS2Scan:
    rt: float | None
    native_id: str
    spectrum_index: int
    scan_number: int | None
    polarity: str
    precursor_mz: float | None = None
    precursor_intensity: float | None = None
    precursor_charge: int | None = None
    precursor_ref: str | None = None
    isolation_width: float | None = None
    precursor_count: int = 0
    precursor_mz_source: str = "missing"
    centroided: bool | None = None
    mz: np.ndarray | None = None
    intensity: np.ndarray | None = None
    intensity_source: str = "unavailable"
    fragment_intensity: float | None = None
    polarity_source: str = "metadata"
    _intensity_polarity: str | None = None


def _scan_number(native_id):
    if native_id.isdecimal():
        return int(native_id)
    match = re.search(r"(?:^|\s)scan=(\d+)(?:\s|$)", native_id)
    return int(match.group(1)) if match else None


def _optional_number(value, *, positive=False):
    if value is None or value == "":
        return None
    number = float(value)
    if not math.isfinite(number) or number < 0 or (positive and number == 0):
        raise ValueError(f"Invalid spectrum metadata number: {value}")
    return number


def _acquisition_range(lower, upper):
    """Bad optional bounds must not discard an otherwise readable spectrum."""
    if lower is None and upper is None:
        return None, False
    try:
        lower, upper = float(lower), float(upper)
        if math.isfinite(lower) and math.isfinite(upper) and 0 <= lower < upper:
            return (lower, upper), False
    except (TypeError, ValueError):
        pass
    return None, True


def _polarity(params):
    if "MS:1000130" in params and "MS:1000129" in params:
        raise ValueError("Conflicting positive and negative polarity metadata")
    return (
        "positive"
        if "MS:1000130" in params
        else "negative"
        if "MS:1000129" in params
        else "unknown"
    )


def _local(element):
    return etree.QName(element).localname if isinstance(element.tag, str) else ""


def _release(element):
    element.clear()
    parent = element.getparent()
    if parent is not None:
        while element.getprevious() is not None:
            del parent[0]


def _params(element, groups):
    params = {}
    for child in element:
        name = _local(child)
        if name == "referenceableParamGroupRef":
            ref = child.get("ref")
            if ref not in groups:
                raise ValueError(f"Unresolved referenceableParamGroup: {ref}")
            params.update(groups[ref])
        elif name == "cvParam":
            params[child.get("accession")] = dict(child.attrib)
    return params


def _centroid(mode, assume_centroid, native_id):
    if mode is None and not assume_centroid:
        raise ValueError(
            f"Unknown centroid status in MS1 spectrum {native_id}; use assume_centroid=True "
            "(--assume-centroid) only if the file is known to be centroided"
        )


def _arrays(mz, intensity, *, profile=False):
    mz = np.ascontiguousarray(mz, dtype=np.float64)
    intensity = np.ascontiguousarray(intensity, dtype=np.float64)
    if mz.ndim != 1 or intensity.ndim != 1 or mz.size != intensity.size:
        raise ValueError("Mismatched m/z and intensity arrays")
    if normalized_arrays(mz, intensity, True, not profile):
        return mz, intensity
    if not np.all(np.isfinite(mz)) or not np.all(np.isfinite(intensity)):
        raise ValueError("Non-finite m/z or intensity")
    if np.any(mz <= 0) or (not profile and np.any(intensity < 0)):
        raise ValueError("m/z must be positive and intensity must be nonnegative")
    if profile:
        # Baseline-corrected profile exports can contain negative samples.
        intensity = np.maximum(intensity, 0)
    positive = intensity > 0
    if not profile and not positive.all():
        mz, intensity = mz[positive], intensity[positive]
    if np.any(mz[1:] < mz[:-1]):
        order = np.argsort(mz, kind="stable")
        mz, intensity = mz[order], intensity[order]
    if profile and np.any(mz[1:] == mz[:-1]):
        starts = np.r_[0, np.flatnonzero(mz[1:] != mz[:-1]) + 1]
        mz, intensity = mz[starts], np.add.reduceat(intensity, starts)
        if not np.all(np.isfinite(intensity)):
            raise ValueError("Intensity overflow while merging duplicate centroids")
    return mz, intensity


def _base64(text):
    return base64.b64decode("".join((text or "").split()), validate=True)


_NUMPRESS = {
    "MS:1002312": ("linear", False),
    "MS:1002313": ("pic", False),
    "MS:1002314": ("slof", False),
    "MS:1002746": ("linear", True),
    "MS:1002747": ("pic", True),
    "MS:1002748": ("slof", True),
}
_COMPRESSION = set(_NUMPRESS) | {"MS:1000574", "MS:1000576"}


def _mzml_binary(element, params):
    compression = set(params) & _COMPRESSION
    unknown = [
        p.get("name", "")
        for k, p in params.items()
        if "compression" in p.get("name", "").lower() and k not in _COMPRESSION
    ]
    if not compression or unknown:
        raise ValueError(f"Unsupported or missing mzML compression: {unknown or list(params)}")
    if "MS:1000576" in compression and len(compression) > 1:
        raise ValueError("Conflicting mzML compression terms")
    numpress = list(set(params) & set(_NUMPRESS))
    if len(numpress) > 1:
        raise ValueError("Conflicting MS-Numpress encodings")
    binary = element.find("{*}binary")
    if binary is None:
        raise ValueError("Missing mzML binary element")
    raw = _base64(binary.text)
    if "MS:1000574" in params or (numpress and _NUMPRESS[numpress[0]][1]):
        raw = zlib.decompress(raw)
    if numpress:
        try:
            import pynumpress
        except ImportError as error:
            raise ValueError(
                "MS-Numpress input requires: pip install 'masscubex[numpress]'"
            ) from error
        decoder = getattr(pynumpress, "decode_" + _NUMPRESS[numpress[0]][0])
        return decoder(np.frombuffer(raw, dtype=np.uint8))
    precisions = set(params) & {"MS:1000521", "MS:1000523"}
    if len(precisions) != 1:
        raise ValueError("mzML arrays require one 32-bit or 64-bit floating-point type")
    # mzML 1.1 specifies little-endian binary arrays.
    return np.frombuffer(raw, dtype="<f4" if "MS:1000521" in params else "<f8")


def _mzml_scan(
    element, groups, assume_centroid, include_ms2=False, decode_ms2=False, index=-1, sum_ms2=False
):
    params = _params(element, groups)
    if "MS:1000511" not in params:
        raise ValueError("mzML spectrum has no MS level")
    level = int(params["MS:1000511"]["value"])
    if level == 2 and include_ms2:
        return _mzml_ms2(element, groups, params, decode_ms2, index, sum_ms2)
    if level != 1:
        return None
    native_id = element.get("id", "")
    mode = False if "MS:1000128" in params else True if "MS:1000127" in params else None
    _centroid(mode, assume_centroid, native_id)
    polarity = _polarity(params)
    scan = element.find("{*}scanList/{*}scan")
    scan_params = _params(scan, groups) if scan is not None else {}
    # A scalar mobility coordinate means frames must be collapsed beforehand.
    if set(scan_params) & {"MS:1002476", "MS:1002815", "MS:1001581"}:
        raise ValueError("Ion-mobility-resolved MS1 scans must be collapsed before extraction")
    rt = _mzml_time(scan_params, native_id, required=True)
    mz, intensity = _mzml_arrays(element, groups, profile=mode is False)
    ranges, invalid = [], False
    # Read each window separately: flattening CV terms loses multiple windows.
    for window in element.findall("{*}scanList/{*}scan/{*}scanWindowList/{*}scanWindow"):
        terms = _params(window, groups)
        bounds, bad = _acquisition_range(
            terms.get("MS:1000501", {}).get("value"),
            terms.get("MS:1000500", {}).get("value"),
        )
        invalid |= bad
        if bounds is not None:
            ranges.append(bounds)
    acquired = (min(r[0] for r in ranges), max(r[1] for r in ranges)) if ranges else None
    return Scan(rt, mz, intensity, polarity, native_id, mode is not False, index, acquired, invalid)


def _mzml_time(scan_params, native_id, required=False):
    time = scan_params.get("MS:1000016")
    if time is None:
        if required:
            raise ValueError(f"MS1 spectrum {native_id} has no scan start time")
        return None
    unit = time.get("unitAccession") or time.get("unitName", "").lower()
    factors = {
        "UO:0000010": 1,
        "UO:0000031": 60,
        "UO:0000032": 3600,
        "UO:0000028": 0.001,
        "second": 1,
        "minute": 60,
        "hour": 3600,
        "millisecond": 0.001,
    }
    if unit not in factors:
        raise ValueError(f"Missing or unsupported scan time unit: {unit}")
    return float(time["value"]) * factors[unit]


def _mzml_arrays(element, groups, profile=False):
    expected = int(element.get("defaultArrayLength", "-1"))
    if expected < 0:
        raise ValueError("Missing or invalid defaultArrayLength")
    arrays = {}
    for array in element.findall("{*}binaryDataArrayList/{*}binaryDataArray"):
        array_params = _params(array, groups)
        keys = set(array_params) & {"MS:1000514", "MS:1000515"}
        if not keys:
            continue
        if len(keys) != 1:
            raise ValueError("Ambiguous binary array type")
        key = keys.pop()
        if key in arrays:
            raise ValueError("Duplicate mzML m/z or intensity array")
        values = _mzml_binary(array, array_params)
        if len(values) != expected:
            raise ValueError(f"Binary array length {len(values)} != defaultArrayLength {expected}")
        arrays[key] = values
    if expected and len(arrays) != 2:
        raise ValueError("Missing m/z or intensity array")
    mz, intensity = _arrays(
        arrays.get("MS:1000514", []), arrays.get("MS:1000515", []), profile=profile
    )
    return mz, intensity


def _ms2_arrays(result, mz, intensity, retain):
    """Score the complete exported spectrum before optional evidence pruning."""
    result.fragment_intensity = float(np.sum(intensity, dtype=np.float64))
    if not math.isfinite(result.fragment_intensity):
        raise ValueError("MS2 summed fragment intensity overflow")
    if retain:
        result.mz, result.intensity = mz, intensity


def _mzml_ms2(element, groups, params, decode, index, sum_ms2=False):
    native_id = element.get("id", "")
    scan = element.find("{*}scanList/{*}scan")
    scan_params = _params(scan, groups) if scan is not None else {}
    mode = False if "MS:1000128" in params else True if "MS:1000127" in params else None
    result = MS2Scan(
        _mzml_time(scan_params, native_id),
        native_id,
        index,
        _scan_number(native_id),
        _polarity(params),
        centroided=mode,
    )
    precursors = element.findall("{*}precursorList/{*}precursor")
    selected = [
        (precursor, ion)
        for precursor in precursors
        for ion in precursor.findall("{*}selectedIonList/{*}selectedIon")
    ]
    result.precursor_count = max(len(precursors), len(selected))
    if result.precursor_count == 1:
        precursor = precursors[0]
        result.precursor_ref = precursor.get("spectrumRef")
        isolation = precursor.find("{*}isolationWindow")
        window = _params(isolation, groups) if isolation is not None else {}
        lower = _optional_number(window.get("MS:1000828", {}).get("value"))
        upper = _optional_number(window.get("MS:1000829", {}).get("value"))
        if lower is not None and upper is not None:
            result.isolation_width = lower + upper
        ions = _params(selected[0][1], groups) if selected else {}
        result.precursor_mz = _optional_number(
            ions.get("MS:1000744", {}).get("value"), positive=True
        )
        result.precursor_mz_source = (
            "selected_ion" if result.precursor_mz is not None else "missing"
        )
        if result.precursor_mz is None and result.isolation_width is not None:
            result.precursor_mz = _optional_number(
                window.get("MS:1000827", {}).get("value"), positive=True
            )
            if result.precursor_mz is not None:
                result.precursor_mz_source = "isolation_target"
        result.precursor_intensity = _optional_number(ions.get("MS:1000042", {}).get("value"))
        charge = ions.get("MS:1000041", {}).get("value")
        result.precursor_charge = int(charge) if charge not in (None, "") else None
    if sum_ms2 or (decode and mode is True):
        mz, intensity = _mzml_arrays(element, groups, profile=mode is False)
        _ms2_arrays(result, mz, intensity, decode)
    return result


def _iter_mzml(stream, assume_centroid, include_ms2=False, decode_ms2=False, sum_ms2=False):
    groups = {}
    index = 0
    root_seen = False
    tags = tuple(
        "{*}" + name
        for name in (
            "mzML",
            "indexedmzML",
            "referenceableParamGroup",
            "spectrum",
            "chromatogram",
            "offset",
            "index",
        )
    )
    events = etree.iterparse(
        stream,
        events=("start", "end"),
        tag=tags,
        resolve_entities=False,
        no_network=True,
        load_dtd=False,
        huge_tree=True,
    )
    for event, element in events:
        tag = _local(element)
        if event == "start":
            if element.getparent() is None:
                if tag not in {"mzML", "indexedmzML"}:
                    raise ValueError("Expected an mzML document")
                root_seen = True
            continue
        if tag == "referenceableParamGroup":
            groups[element.get("id")] = _params(element, groups)
            _release(element)
        elif tag == "spectrum":
            scan = _mzml_scan(
                element, groups, assume_centroid, include_ms2, decode_ms2, index, sum_ms2
            )
            index += 1
            _release(element)
            if scan is not None:
                yield scan
        elif tag in {"chromatogram", "offset", "index"}:
            _release(element)
    if not root_seen:
        raise ValueError("Expected an mzML document")


_DURATION = re.compile(
    r"P(?:(\d+(?:\.\d+)?)D)?(?:T(?:(\d+(?:\.\d+)?)H)?(?:(\d+(?:\.\d+)?)M)?(?:(\d+(?:\.\d+)?)S)?)?"
)


def _duration(value):
    match = _DURATION.fullmatch(value)
    if match is None or not any(match.groups()) or value.endswith("T"):
        raise ValueError(f"Invalid mzXML retentionTime: {value}")
    days, hours, minutes, seconds = (float(v or 0) for v in match.groups())
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def _mzxml_scan(
    attrs,
    peaks,
    global_centroid,
    assume_centroid,
    include_ms2=False,
    decode_ms2=False,
    index=-1,
    precursors=(),
    sum_ms2=False,
):
    level = int(attrs["msLevel"])
    if level == 2 and include_ms2:
        return _mzxml_ms2(attrs, peaks, precursors, global_centroid, decode_ms2, index, sum_ms2)
    if level != 1:
        return None
    native_id = attrs.get("num", "")
    flag = attrs.get("centroided", global_centroid)
    if flag not in (None, "0", "1", "true", "false"):
        raise ValueError(f"Invalid centroided value: {flag}")
    _centroid(None if flag is None else flag in {"1", "true"}, assume_centroid, native_id)
    rt = _duration(attrs.get("retentionTime", ""))
    profile = flag in {"0", "false"}
    mz, intensity = _mzxml_arrays(attrs, peaks, profile)
    # startMz/endMz describe the acquisition; lowMz/highMz are observed extrema.
    acquired, invalid = _acquisition_range(attrs.get("startMz"), attrs.get("endMz"))
    return Scan(
        rt,
        mz,
        intensity,
        {"+": "positive", "-": "negative"}.get(attrs.get("polarity"), "unknown"),
        native_id,
        not profile,
        index,
        acquired,
        invalid,
    )


def _mzxml_arrays(attrs, peaks, profile=False):
    count = int(attrs.get("peaksCount", "-1"))
    if count < 0:
        raise ValueError("Missing or invalid peaksCount")
    if peaks is None or (count == 0 and not (peaks.text or "").strip()):
        if count:
            raise ValueError("Missing mzXML peaks")
        values = np.empty((0, 2))
    else:
        compression = peaks.get("compressionType", "none")
        if compression not in {"none", "zlib"}:
            raise ValueError(f"Unsupported mzXML compression: {compression}")
        order = peaks.get("byteOrder", "network")
        if order not in {"network", "big", "little"}:
            raise ValueError(f"Unsupported mzXML byteOrder: {order}")
        precision = peaks.get("precision")
        if precision not in {"32", "64"}:
            raise ValueError(f"Unsupported mzXML precision: {precision}")
        content = peaks.get("pairOrder", peaks.get("contentType", "m/z-int"))
        if content not in {"m/z-int", "m/z-intensity"}:
            raise ValueError(f"Unsupported mzXML peak ordering: {content}")
        raw = _base64(peaks.text)
        if compression == "zlib":
            raw = zlib.decompress(raw)
        dtype = ("<" if order == "little" else ">") + ("f4" if precision == "32" else "f8")
        values = np.frombuffer(raw, dtype=dtype)
        if values.size != 2 * count:
            raise ValueError("mzXML binary array length does not match peaksCount")
        values = values.reshape(-1, 2)
    mz, intensity = _arrays(values[:, 0], values[:, 1], profile=profile)
    return mz, intensity


def _mzxml_ms2(attrs, peaks, precursors, global_centroid, decode, index, sum_ms2=False):
    flag = attrs.get("centroided", global_centroid)
    if flag not in (None, "0", "1", "true", "false"):
        raise ValueError(f"Invalid centroided value: {flag}")
    mode = None if flag is None else flag in {"1", "true"}
    native_id = attrs.get("num", "")
    time = attrs.get("retentionTime")
    result = MS2Scan(
        _duration(time) if time else None,
        native_id,
        index,
        _scan_number(native_id),
        {"+": "positive", "-": "negative"}.get(attrs.get("polarity"), "unknown"),
        precursor_count=len(precursors),
        centroided=mode,
    )
    if len(precursors) == 1:
        precursor, mass = precursors[0]
        result.precursor_mz = _optional_number(mass, positive=True)
        result.precursor_mz_source = (
            "selected_ion" if result.precursor_mz is not None else "missing"
        )
        result.precursor_intensity = _optional_number(precursor.get("precursorIntensity"))
        result.precursor_ref = precursor.get("precursorScanNum")
        result.isolation_width = _optional_number(precursor.get("windowWideness"))
        charge = precursor.get("precursorCharge")
        result.precursor_charge = int(charge) if charge not in (None, "") else None
    if sum_ms2 or (decode and mode is True):
        mz, intensity = _mzxml_arrays(attrs, peaks, profile=mode is False)
        _ms2_arrays(result, mz, intensity, decode)
    return result


def _iter_mzxml(stream, assume_centroid, include_ms2=False, decode_ms2=False, sum_ms2=False):
    stack = []
    index = 0
    global_centroid = None
    root_seen = False
    tags = tuple(
        "{*}" + name
        for name in (
            "mzXML",
            "msRun",
            "scan",
            "peaks",
            "precursorMz",
            "dataProcessing",
            "offset",
            "index",
        )
    )
    events = etree.iterparse(
        stream,
        events=("start", "end"),
        tag=tags,
        resolve_entities=False,
        no_network=True,
        load_dtd=False,
        huge_tree=True,
    )
    for event, element in events:
        tag = _local(element)
        if event == "start":
            if element.getparent() is None:
                if tag not in {"mzXML", "msRun"}:
                    raise ValueError("Expected an mzXML document")
                root_seen = True
            if tag == "dataProcessing" and element.get("centroided") is not None:
                global_centroid = element.get("centroided")
            elif tag == "scan":
                if "msLevel" not in element.attrib:
                    raise ValueError("mzXML scan has no MS level")
                stack.append([dict(element.attrib), False, index, []])
                index += 1
        elif tag == "precursorMz":
            if stack:
                stack[-1][3].append((dict(element.attrib), element.text))
        elif tag == "peaks":
            if not stack:
                raise ValueError("mzXML peaks outside a scan")
            if stack[-1][1]:
                raise ValueError("Duplicate mzXML peaks element")
            scan = _mzxml_scan(
                stack[-1][0],
                element,
                global_centroid,
                assume_centroid,
                include_ms2,
                decode_ms2,
                stack[-1][2],
                stack[-1][3],
                sum_ms2,
            )
            stack[-1][1] = True
            _release(element)
            # Emit on peaks end, before nested MS2 scans can change document order.
            if scan is not None:
                yield scan
        elif tag == "scan":
            attrs, emitted, scan_index, precursors = stack.pop()
            if not emitted:
                scan = _mzxml_scan(
                    attrs,
                    None,
                    global_centroid,
                    assume_centroid,
                    include_ms2,
                    decode_ms2,
                    scan_index,
                    precursors,
                    sum_ms2,
                )
                if scan is not None:
                    yield scan
            _release(element)
        elif tag in {"offset", "index"}:
            _release(element)
    if not root_seen:
        raise ValueError("Expected an mzXML document")


def iter_ms1(path, *, assume_centroid=False):
    """Yield normalized MS1 scans in acquisition order, times in seconds.

    Explicit profile spectra retain their samples for mass detection.
    ``assume_centroid`` only allows missing centroid metadata; it never
    re-labels explicit profile data. Non-MS1 spectra skip binary decoding.
    """
    yield from iter_scans(path, assume_centroid=assume_centroid, include_ms2=False)


def iter_scans(path, *, assume_centroid=False, include_ms2=True, decode_ms2=False, sum_ms2=False):
    """Yield MS1 arrays and MS2 metadata in acquisition order.

    MS2 arrays remain undecoded unless requested for fragment evidence or
    summed-intensity ranking. Summing alone immediately discards the arrays.
    Native scan numbers and zero-based spectrum indices are separate identifiers.
    """
    path = Path(path)
    kind = file_format(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    reader = _iter_mzml if kind == "mzml" else _iter_mzxml
    try:
        with opener(path, "rb") as stream:
            for scan in reader(stream, assume_centroid, include_ms2, decode_ms2, sum_ms2):
                if scan.rt is not None and (not math.isfinite(scan.rt) or scan.rt < 0):
                    raise ValueError(f"Invalid retention time in scan {scan.native_id}")
                yield scan
    except (ValueError, KeyError, etree.XMLSyntaxError, zlib.error) as error:
        raise ValueError(f"{path}: {error}") from error
