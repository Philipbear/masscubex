"""Inspect acquisition metadata without decoding or retaining binary spectra.

DIA evidence comes from declared acquisition modes or repeated, densely tiled
isolation windows. Neither a filename nor one wide isolation window is proof.
Missing evidence remains unknown rather than being labelled DDA.
"""

from collections import Counter, deque
from dataclasses import asdict, dataclass, field
import math
from pathlib import Path
import re
import unicodedata

from lxml import etree

from .formats import file_format, open_xml
from .spectrum_types import NON_MS_SPECTRUM_TERMS


# PSI-MS 4.2.14. Acquisition aspects 3216–3218 describe full-range, scanning
# quadrupole and sequential-range dissociation; ion mobility alone is not DIA.
DIA_TERMS = {
    "MS:1003215": "dia",
    "MS:1003216": "all_ion",
    "MS:1003217": "scanning_window",
    "MS:1003218": "windowed",
    "MS:1003224": "windowed",
    "MS:1003225": "dia_pasef",
    "MS:1003226": "all_ion_mobility",
    "MS:1003227": "all_ion",
    "MS:1003228": "scanning_window",
}
DDA_TERM = "MS:1003221"
_MODE = re.compile(
    r"\b(?:data[ -]?independent(?:\s+acquisition)?|dia(?:pasef)?|swath|"
    r"all[ -]?ions?\s+fragmentation|aif|(?:hd)?ms(?:\^|\()?e\)?|bbcid|sonar)\b",
    re.I,
)
_MODE_FIELD = re.compile(r"acquisition|scan.?type|experiment.?type|method.?type|^dia$", re.I)
_NEGATED = re.compile(r"\b(?:no|not|non)[ -]+(?:dia|swath|aif|data[ -]?independent)\b", re.I)
_ALL_ION_FILTER = re.compile(r"\b(?:full\s+)?ms2\s+0+(?:\.0+)?\s*@", re.I)
_FALSE = {"false", "no", "off", "0"}
_VENDORS = (
    ("Thermo", re.compile(r"thermo|finnigan|orbitrap|exactive|\bltq\b", re.I)),
    ("SCIEX", re.compile(r"sciex|applied biosystems|tripletof|qstar|qtrap", re.I)),
    ("Agilent", re.compile(r"agilent|\b65\d\d\b.*(?:tof|q-tof)", re.I)),
    ("Waters", re.compile(r"waters|synapt|xevo|q-tof ultima", re.I)),
    ("Bruker", re.compile(r"bruker|micr[o0]tof|maxis|impact|tims", re.I)),
    ("Shimadzu", re.compile(r"shimadzu|lcms[- ]?(?:90|80)|gcms[- ]?qp", re.I)),
)


@dataclass(frozen=True)
class AcquisitionInfo:
    mode: str
    confidence: str
    evidence: tuple[str, ...]
    dia_families: tuple[str, ...]
    vendor: str
    instrument_models: tuple[str, ...]
    spectra_inspected: int
    ms_levels: dict
    inspection_complete: bool
    binary_arrays_decoded: bool = False
    detector_version: int = 2
    chromatograms_inspected: int = 0
    non_ms_spectra_inspected: int = 0

    @property
    def is_dia(self):
        return self.mode in {"dia", "mixed_dia"}

    def to_dict(self):
        result = asdict(self)
        for name in ("evidence", "dia_families", "instrument_models"):
            result[name] = list(result[name])
        return result


@dataclass
class _Spectrum:
    level: int | str = 0
    polarity: str = "unknown"
    native_id: str = ""
    targets: list = field(default_factory=list)
    selected_ions: list = field(default_factory=list)
    precursor_count: int = 0
    lower: float | None = None
    upper: float | None = None
    width: float | None = None
    scan_lower: float | None = None
    scan_upper: float | None = None
    scan_ranges: list = field(default_factory=list)
    emitted: bool = False
    non_ms: bool = False


def _number(value):
    """Optional detector metadata must not be coerced into a false window."""
    try:
        value = float(value)
    except (ValueError, TypeError):
        return None
    return value if math.isfinite(value) else None


def _tiled(windows, survey=None):
    """Require several different windows that densely span a substantial range."""
    unique = sorted(set(windows))
    if len(unique) < 2:
        return False
    first, last = unique[0][0], max(b for _, b in unique)
    if last - first < 100:
        return False
    covered, end = 0.0, first
    for left, right in unique:
        covered += max(0.0, right - max(left, end))
        end = max(end, right)
    if covered / (last - first) < 0.90:
        return False
    if len(unique) >= 6:
        return True
    # A small number of broad tiles needs additional evidence: together they
    # must cover the acquired MS1 range, rather than a selected PRM interval.
    return survey is not None and (
        max(0.0, min(last, survey[1]) - max(first, survey[0])) / (survey[1] - survey[0]) >= 0.90
    )


class _WindowCycles:
    """Bounded rolling evidence, separate for each polarity.

    Compare complete MS1-delimited cycles. A second rolling check covers
    MS2-only exports, where no survey scan supplies cycle boundaries.
    """

    def __init__(self):
        self.current = []
        self.invalid = False
        self.cycles = deque(maxlen=3)
        self.recent = deque(maxlen=384)
        self.survey = None

    def boundary(self):
        if self.invalid:
            self.cycles.clear()
        elif self.current:
            self.cycles.append(tuple(self.current))
        self.current.clear()
        self.invalid = False
        if len(self.cycles) == 3:
            a, b, c = self.cycles
            return a == b == c and _tiled(a, self.survey)
        return False

    def add(self, window):
        if not self.invalid and len(self.current) < 256:
            self.current.append(window)
        elif not self.invalid:
            # An oversized cycle is unsupported evidence, not a truncated fit.
            self.current.clear()
            self.cycles.clear()
            self.invalid = True
        self.recent.append(window)
        count = len(self.recent)
        if count < 6 or count % 3:
            return False
        recent = list(self.recent)
        # Repeated blocks must be identical, cover >=100 Da and densely tile
        # that span. Repeated sparse PRM targets therefore do not qualify.
        for size in range(2, min(128, count // 3) + 1):
            a, b, c = recent[-3 * size : -2 * size], recent[-2 * size : -size], recent[-size:]
            if a == b == c and _tiled(a, self.survey):
                return True
        return False

    def interrupt(self):
        """Do not infer cycles by silently dropping intervening MS2 scans."""
        self.invalid = True
        self.current.clear()
        self.cycles.clear()
        self.recent.clear()


class _DetectedDIA(Exception):
    pass


class _MetadataTarget:
    """lxml target parser: binary character data are discarded immediately."""

    def __init__(self, kind):
        self.kind = kind
        self.path = []
        self.spectra = []
        self.groups = {}
        self.group = None
        self.precursor = None
        self.precursor_text = []
        self.levels = Counter()
        self.models = []
        self.manufacturers = []
        self.evidence = []
        self.families = set()
        self.dda = False
        self.cycles = {}
        self.ms1_ranges = {}
        self.full_range_windows = {}
        self.spectra_inspected = 0
        self.chromatograms_inspected = 0
        self.non_ms_spectra_inspected = 0
        self.root_seen = False

    def _dia(self, family, evidence):
        self.families.add(family)
        if evidence not in self.evidence:
            self.evidence.append(evidence)
        # Finish file-level headers to retain instrument identity. Spectrum
        # evidence can stop immediately, before expensive binary text is read.
        if self.spectra:
            raise _DetectedDIA

    def _mode(self, name, value, accession=None):
        if str(value or "").strip().lower() in _FALSE:
            return
        if accession in DIA_TERMS:
            self._dia(DIA_TERMS[accession], f"PSI-MS {accession}: {name}")
        elif accession == DDA_TERM:
            self.dda = True
        elif _MODE_FIELD.search(name):
            # Structured boolean flags can put the mode in the field name.
            is_flag = str(value).strip().lower() in {"true", "yes", "on", "1"}
            text = name if is_flag else value or name
            text = unicodedata.normalize("NFKC", text)
            if _NEGATED.search(text):
                return
            if _MODE.search(text):
                self._dia("declared_dia", f"Declared {name}: {value}"[:300])
            elif re.search(r"\bdda\b|data[ -]?dependent", text, re.I):
                self.dda = True

    def _param(self, attrs, parent):
        accession, name, value = (
            attrs.get("accession"),
            attrs.get("name", ""),
            attrs.get("value", ""),
        )
        if self.group is not None:
            self.group[1].append(dict(attrs))
            return
        if self.spectra and parent == "spectrum" and accession in NON_MS_SPECTRUM_TERMS:
            self.spectra[-1].non_ms = True
        if self.spectra and self.spectra[-1].non_ms:
            return
        # Software names, filenames and free-text study comments do not label
        # the acquisition. Only acquisition terms/fields are authoritative.
        if parent not in {"software", "sourceFile", "contact", "sample"}:
            self._mode(name, value, accession)
        if parent == "instrumentConfiguration":
            model = value if accession == "MS:1000031" else name
            if model and not re.search(r"serial number|instrument attribute", name, re.I):
                if model not in self.models:
                    self.models.append(model)
        elif name.lower() in {"msmanufacturer", "manufacturer"} and not self.spectra:
            self.manufacturers.append(value)
        if not self.spectra:
            return
        scan = self.spectra[-1]
        if accession == "MS:1000511":
            scan.level = value
        elif accession in {"MS:1000130", "MS:1000129"}:
            scan.polarity = "positive" if accession == "MS:1000130" else "negative"
        elif accession == "MS:1000827":
            scan.targets.append(_number(value))
        elif accession == "MS:1000744" and parent == "selectedIon":
            scan.selected_ions.append(_number(value))
        elif accession == "MS:1000828":
            scan.lower = _number(value)
        elif accession == "MS:1000829":
            scan.upper = _number(value)
        elif accession == "MS:1000501":
            scan.scan_lower = _number(value)
        elif accession == "MS:1000500":
            scan.scan_upper = _number(value)
        elif accession == "MS:1000512" and _ALL_ION_FILTER.search(value):
            self._dia("all_ion", f"All-ion scan filter: {value}"[:300])

    def start(self, tag, attrs):
        tag = tag.rsplit("}", 1)[-1]
        parent = self.path[-1] if self.path else None
        if not self.path:
            valid = {"mzML", "indexedmzML"} if self.kind == "mzml" else {"mzXML", "msRun"}
            if tag not in valid:
                raise ValueError(f"Expected {self.kind} XML, found {tag}")
            self.root_seen = True
        self.path.append(tag)
        if tag == "referenceableParamGroup":
            self.group = (attrs.get("id"), [])
        elif tag == "referenceableParamGroupRef":
            reference = attrs.get("ref")
            if reference not in self.groups:
                raise ValueError(f"Unresolved referenceableParamGroup: {reference}")
            for param in self.groups[reference]:
                self._param(param, parent)
        elif tag in {"cvParam", "userParam"}:
            self._param(attrs, parent)
        elif tag == "msManufacturer":
            self.manufacturers.append(attrs.get("value", ""))
        elif tag == "msModel":
            self.models.append(attrs.get("value", ""))
        elif tag == "chromatogram":
            self.chromatograms_inspected += 1
        elif tag == "spectrum" or (self.kind == "mzxml" and tag == "scan"):
            if self.families:
                raise _DetectedDIA
            if self.kind == "mzxml" and self.spectra:
                self._emit(self.spectra[-1])  # Parent scan precedes nested MS2.
            scan = _Spectrum(native_id=attrs.get("id", attrs.get("num", "")))
            if self.kind == "mzxml":
                scan.level = int(attrs.get("msLevel", 0))
                scan.polarity = {"+": "positive", "-": "negative"}.get(
                    attrs.get("polarity"), "unknown"
                )
                scan.scan_lower = _number(attrs.get("startMz"))
                scan.scan_upper = _number(attrs.get("endMz"))
            self.spectra.append(scan)
            if self.kind == "mzxml":
                self._mode("scan type", attrs.get("scanType", ""))
                scan_filter = attrs.get("filterLine", "")
                if _ALL_ION_FILTER.search(scan_filter):
                    self._dia("all_ion", f"All-ion scan filter: {scan_filter}"[:300])
        elif self.kind == "mzxml" and tag == "precursorMz" and self.spectra:
            self.spectra[-1].precursor_count += 1
            self.precursor = dict(attrs)
            self.precursor_text = []
        elif tag == "precursor" and self.spectra:
            self.spectra[-1].precursor_count += 1

    def data(self, text):
        if self.precursor is not None and self.path[-1] == "precursorMz":
            if sum(map(len, self.precursor_text)) + len(text) > 128:
                raise ValueError("Invalid mzXML precursor mass text")
            self.precursor_text.append(text)
        # In particular, do not retain <binary> or <peaks> payloads.

    def _emit(self, scan):
        if scan.emitted:
            return
        scan.emitted = True
        self.spectra_inspected += 1
        if scan.non_ms:
            self.non_ms_spectra_inspected += 1
            return
        scan.level = int(scan.level)
        self.levels[scan.level] += 1
        cycle = self.cycles.setdefault(scan.polarity, _WindowCycles())
        if scan.level == 1:
            if scan.scan_ranges:
                scan.scan_lower = min(a for a, _ in scan.scan_ranges)
                scan.scan_upper = max(b for _, b in scan.scan_ranges)
            if (
                scan.scan_lower is not None
                and scan.scan_upper is not None
                and scan.scan_upper - scan.scan_lower >= 100
            ):
                # A narrow lock-mass survey must not redefine the full range.
                self.ms1_ranges[scan.polarity] = (scan.scan_lower, scan.scan_upper)
                cycle.survey = self.ms1_ranges[scan.polarity]
            if cycle.boundary():
                self._dia("windowed", "Three identical MS1-delimited cycles densely tile >=100 Da")
        elif scan.level == 2:
            targets = scan.targets or scan.selected_ions
            valid = [x for x in targets if x is not None and x > 0]
            width = scan.width
            if width is None and scan.lower is not None and scan.upper is not None:
                width = scan.lower + scan.upper
            offsets_valid = all(x is None or x >= 0 for x in (scan.lower, scan.upper))
            if (
                len(valid) == 1
                and scan.precursor_count <= 1
                and width is not None
                and width > 0
                and offsets_valid
            ):
                lower = valid[0] - (scan.lower if scan.lower is not None else width / 2)
                upper = valid[0] + (scan.upper if scan.upper is not None else width / 2)
                if lower >= 0 and upper > lower:
                    window = (round(lower, 4), round(upper, 4))
                    survey = self.ms1_ranges.get(scan.polarity)
                    if survey is not None:
                        overlap = max(0.0, min(upper, survey[1]) - max(lower, survey[0]))
                        if overlap / (survey[1] - survey[0]) >= 0.90:
                            previous, count = self.full_range_windows.get(scan.polarity, (None, 0))
                            count = count + 1 if window == previous else 1
                            self.full_range_windows[scan.polarity] = (window, count)
                            if count >= 3:
                                self._dia(
                                    "all_ion",
                                    f"Three identical MS2 isolation windows {window} cover >=90% of "
                                    f"the acquired MS1 range {survey}",
                                )
                        else:
                            self.full_range_windows.pop(scan.polarity, None)
                    if cycle.add(window):
                        self._dia(
                            "windowed",
                            "Three repeated isolation-window blocks densely tile >=100 Da",
                        )
                    return
            cycle.interrupt()
            self.full_range_windows.pop(scan.polarity, None)

    def end(self, tag):
        tag = tag.rsplit("}", 1)[-1]
        if tag == "referenceableParamGroup":
            self.groups[self.group[0]] = self.group[1]
            self.group = None
        elif self.kind == "mzxml" and tag == "precursorMz" and self.precursor is not None:
            scan = self.spectra[-1]
            scan.targets.append(_number("".join(self.precursor_text)))
            scan.width = _number(self.precursor.get("windowWideness"))
            self.precursor = None
            self.precursor_text.clear()
        elif tag == "scanWindow" and self.spectra:
            scan = self.spectra[-1]
            if scan.scan_lower is not None and scan.scan_upper is not None:
                if 0 <= scan.scan_lower < scan.scan_upper:
                    scan.scan_ranges.append((scan.scan_lower, scan.scan_upper))
            scan.scan_lower = scan.scan_upper = None
        elif tag == "spectrum" or (self.kind == "mzxml" and tag == "scan"):
            self._emit(self.spectra[-1])
            self.spectra.pop()
        self.path.pop()

    def close(self):
        return None

    def result(self, complete):
        identity = " ".join(self.manufacturers + self.models)
        vendor = next((name for name, pattern in _VENDORS if pattern.search(identity)), "unknown")
        if self.families:
            mode, confidence = ("mixed_dia" if self.dda else "dia"), "high"
        elif self.dda:
            mode, confidence = "dda", "declared"
        elif self.non_ms_spectra_inspected and not self.levels:
            mode, confidence = "non_ms_only", "observed"
        elif not self.spectra_inspected and self.chromatograms_inspected:
            mode, confidence = "chromatogram_only", "observed"
        elif self.levels and set(self.levels) == {1}:
            mode, confidence = "ms1_only", "observed"
        else:
            mode, confidence = "unknown", "insufficient_metadata"
        return AcquisitionInfo(
            mode,
            confidence,
            tuple(self.evidence),
            tuple(sorted(self.families)),
            vendor,
            tuple(dict.fromkeys(self.models)),
            self.spectra_inspected,
            {str(k): v for k, v in sorted(self.levels.items())},
            complete,
            chromatograms_inspected=self.chromatograms_inspected,
            non_ms_spectra_inspected=self.non_ms_spectra_inspected,
        )


def inspect_acquisition(path):
    """Read all metadata for non-DIA files; stop early once DIA is established.

    A target parser discards binary text as it streams. The complete pass for
    other files also detects acquisition switches late in a run. No spectra or
    binary arrays are retained, and no internet/CV download is needed at runtime.
    """
    path = Path(path)
    target = _MetadataTarget(file_format(path))
    parser = etree.XMLParser(
        target=target,
        resolve_entities=False,
        no_network=True,
        load_dtd=False,
        huge_tree=True,
    )
    try:
        with open_xml(path) as stream:
            while chunk := stream.read(1024 * 1024):
                parser.feed(chunk)
            parser.close()
    except _DetectedDIA:
        return target.result(False)
    except (etree.XMLSyntaxError, ValueError) as error:
        raise ValueError(f"Acquisition metadata: {path}: {error}") from error
    if not target.root_seen:
        raise ValueError(f"Acquisition metadata: empty XML document: {path}")
    # Header-only documents can declare DIA without any spectra.
    return target.result(True)
