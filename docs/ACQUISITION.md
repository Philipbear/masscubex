# Acquisition screening

MassCubeX excludes confidently recognized DIA files before calibration, binary
spectrum decoding, MS2 seeding and feature-table creation. This is an input
acceptance policy; the inherited MassCube feature algorithms are unchanged.

```bash
masscubex extract /data/study -o features --resume
masscubex extract sample.mzML -o features --dia-policy error
# Explicit override for a separate workflow that accepts DIA:
masscubex extract sample.mzML -o features --dia-policy process
```

The default is `--dia-policy skip`. A DIA exclusion is a successful policy
decision, not a corrupt-file error. Files that contain any recognized DIA
segment are excluded, including mixed DDA/DIA acquisitions. No DIA deconvolution
or assignment of broad-window fragments to individual MS1 features is performed.

## Evidence

The detector accepts these sources of evidence:

1. **Acquisition declarations.** PSI-MS controlled terms for DIA and its
   full-range, sequential-window and scanning-quadrupole variants; or structured
   acquisition/scan-type fields declaring DIA, SWATH, AIF, MSᴱ, HDMSᴱ, bbCID,
   SONAR or diaPASEF. False/off declarations and negated mode names are ignored.
   File names, sample descriptions and software names are not acquisition evidence.
   Referenceable parameter groups only apply when referenced.
2. **All-ion scan filters.** A Thermo MS2 scan filter with a zero isolation
   target, such as `Full ms2 0.00@hcd...`, declares all-ion fragmentation.
3. **Repeated isolation-window cycles.** Three identical ordered blocks must
   contain distinct windows whose union covers at least 90% of a span of at
   least 100 Da. Six distinct windows suffice. Two to five windows additionally
   need to cover at least 90% of the acquired MS1 scan range. Both MS1-delimited
   cycles and a rolling sequence are checked; the latter supports MS2-only
   exports and repeating combinations of staggered windows. Narrow and
   overlapping windows are supported. Sparse repeated PRM targets do not qualify.
4. **Repeated full-range fragmentation.** Three identical MS2 isolation windows
   each covering at least 90% of an acquired MS1 range of at least 100 Da are
   classified as all-ion acquisition. This recognizes exports that encode
   all-ion fragmentation as a precursor centered on the full acquisition range.
   A narrow lock-mass MS1 scan does not replace that range.

An isolation window uses its target and lower/upper offsets; mzXML uses
`precursorMz` and `windowWideness`. A single selected-ion mass can supply a
missing mzML isolation target when offsets are available. Negative/nonfinite
offsets, missing widths and multiple precursors do not supply structural evidence.
An intervening MS2 with unusable window metadata breaks the pattern; it is not
silently dropped. Evidence is tracked separately for each polarity. There is
**no single isolation-width cutoff**, including no 4 Da cutoff.

The vocabulary mapping uses [PSI-MS 4.2.14](https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo):
`MS:1003215`–`MS:1003218` and `MS:1003224`–`MS:1003228`. Ion mobility alone
(`MS:1003219`) is not DIA evidence. Explicit DDA is `MS:1003221`.
The repeated-window rules above are MassCubeX heuristics, not PSI-MS requirements.
For the underlying all-ion principle, Waters describes
[alternating low/high energy acquisition without precursor selection](https://support.waters.com/KB_Inst/Mass_Spectrometry/WKB84530_What_is_MSe_Acquisition_Method_on_QTof).

## Streaming and uncertainty

A target XML parser discards binary character data as it streams; it never
base64-decodes spectra or builds whole-spectrum XML trees. Pattern storage is
bounded: three cycles of at most 256 windows and a 384-window rolling history
per polarity. The rolling check considers blocks up to 128 windows. Oversized
cycles are invalidated, not fitted after truncation. Compression is recognized
from gzip magic bytes, including files with an incorrect `.gz` suffix.

For files without conclusive DIA evidence, screening reads the complete XML
document so that a late acquisition switch can be found. For confirmed DIA it
stops early; `inspection_complete: false` records that the remainder, including
binary integrity, was not validated. Screening time is recorded separately as
`acquisition_inspection` in the successful file's timing sidecar.

Modes are `dia`, `mixed_dia`, explicitly declared `dda`, observed `ms1_only`,
`chromatogram_only`, or `unknown`. **Unknown is not proof of DDA.** Unknown files
are processed with their uncertainty retained in reports. Exports that omit
mode declarations and usable isolation metadata can conceal DIA. Irregular
adaptive DIA, very long cycles, multiplexed precursors, and all-ion exports
with neither an all-ion declaration nor a full-range isolation window may
remain unknown. A converted MS1-only file cannot establish the acquisition mode
of omitted MS2 spectra. Vendor identity alone never determines acquisition mode.
Dense targeted window schedules can also resemble DIA; mode declarations remain
stronger evidence than inferred window patterns.

Chromatogram-only MRM/SRM exports have no mass spectra and receive a specific
file error before calibration. Mobility-resolved MS1 remains unsupported until
frames are collapsed explicitly. Invalid XML, binary encodings and arrays remain
file errors; no partial feature table is published.

Agilent MassHunter `scanId=N` is retained as its native numeric scan ID. Compound
identifiers such as SCIEX cycle/experiment tuples are not converted into invented
scan numbers; JSON identifier fallbacks preserve their native IDs and source
spectrum indices, as described in [output columns](USAGE.md#output-columns).

## Reports and restart behavior

Each CLI run creates:

| Report | Contents |
| --- | --- |
| `skipped.tsv` | Excluded inputs, reason, acquisition mode, vendor and evidence |
| `errors.tsv` | Corrupt/unsupported/missing inputs and error messages |
| `files.jsonl` | Every outcome, including acquisition details for successful and DIA-excluded inputs |
| `summary.json` | Totals, acquisition-mode counts and skip-reason counts |

`dia_acquisition` and `already_complete` are separate skip reasons. Resume hits
stay in JSONL and summary counts; they are omitted from `skipped.tsv`. A DIA
skip has no new feature table and `output: null`. Existing tables are preserved
if an input is now excluded; use reports to select accepted tables downstream.
Code/configuration fingerprints invalidate older completion markers, so old
results cannot bypass the new screening policy under `--resume`.

The Python `process_file` API returns the same acquisition/skip information;
the CLI persists it in batch reports. File errors do not stop later files, and
the CLI returns a failure exit code after completing a batch with file errors.

