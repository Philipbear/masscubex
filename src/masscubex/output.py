"""Bounded feature-table writers; round measurements only at serialization."""

import csv
from dataclasses import fields
from functools import lru_cache
from itertools import islice
import json
import math
import os
from pathlib import Path
from typing import get_args

from .features import Feature


ALL_FIELDS = fields(Feature)
FEATURE_FIELDS = [f for f in ALL_FIELDS if f.metadata.get("tsv", True)]
FEATURE_COLUMNS = ["feature_id"] + [f.name for f in FEATURE_FIELDS]
FLOAT_COLUMNS = {f.name for f in FEATURE_FIELDS if float in (get_args(f.type) or (f.type,))}
# Historical TSVs must still parse as typed values, including removed columns.
LEGACY_FLOAT_COLUMNS = {f.name for f in ALL_FIELDS if float in (get_args(f.type) or (f.type,))}
NESTED_COLUMNS = {
    "ms2_scan_numbers",
    "ms2_native_ids",
    "ms2_spectrum_indices",
    "ion_relations",
}
BOOLEAN_COLUMNS = {f.name for f in FEATURE_FIELDS if bool in (get_args(f.type) or (f.type,))}
ROW_GROUP_SIZE = 8192
TABLE_SCHEMA_VERSION = 10
# Named policies avoid accidentally rounding an identifier or a new measurement
# as though it were an intensity. Both writers use Python binary64 half-even
# rounding, including halfway cases; the in-memory Feature is never modified.
DECIMAL_PLACES = {
    "mz": 5,
    "rt_apex_seconds": 2,
    "rt_start_seconds": 2,
    "rt_end_seconds": 2,
    "peak_height": 1,
    "peak_area": 1,
    "apex_noise_threshold": 1,
    "gaussian_similarity": 4,
}
RELATION_DECIMAL_PLACES = {
    "neutral_mass": 5,
    "mass_error_ppm": 2,
    "correlation": 4,
    "fragment_mz": 5,
    "fragment_relative_intensity": 1,
}
OUTPUT_PRECISION = {
    "columns": DECIMAL_PLACES,
    "ion_relations": RELATION_DECIMAL_PLACES,
    "retention_time_unit": "second",
    "rounding": "half_even",
}
RELATION_COLUMNS = frozenset(
    {
        "relation",
        "status",
        "partner_feature_id",
        "ion",
        "partner_ion",
        "neutral_mass",
        "mass_error_ppm",
        "correlation",
        "shared_scans",
        "parent_feature_id",
        "fragment_feature_id",
        "parent_ion",
        "fragment_ion",
        "role",
        "ms2_supported",
        "evidence_ms2_scan_number",
        "evidence_ms2_native_id",
        "evidence_ms2_spectrum_index",
        "fragment_mz",
        "fragment_relative_intensity",
        "isotope_state",
    }
)


def resolve_output(source, output_path, output_format):
    if output_format is not None and output_format not in {"parquet", "tsv"}:
        raise ValueError("Output format must be parquet or tsv")
    if output_path is None:
        output_format = output_format or "parquet"
        output_path = str(source) + ".features." + output_format
    path = Path(output_path).resolve()
    inferred = path.suffix.lower().lstrip(".")
    if inferred not in {"parquet", "tsv"}:
        raise ValueError("Output path must end with .parquet or .tsv")
    if output_format is not None and inferred != output_format:
        raise ValueError("Output format does not match the output path extension")
    return path, inferred


@lru_cache(maxsize=1)
def parquet_schema():
    import pyarrow as pa

    # Every possible evidence field has a fixed type, including empty files and
    # row groups with only one relationship kind. Absent evidence is null.
    relation_fields = {
        "relation": pa.string(),
        "status": pa.string(),
        "partner_feature_id": pa.int64(),
        "ion": pa.string(),
        "partner_ion": pa.string(),
        "neutral_mass": pa.float64(),
        "mass_error_ppm": pa.float64(),
        "correlation": pa.float64(),
        "shared_scans": pa.int64(),
        "parent_feature_id": pa.int64(),
        "fragment_feature_id": pa.int64(),
        "parent_ion": pa.string(),
        "fragment_ion": pa.string(),
        "role": pa.string(),
        "ms2_supported": pa.bool_(),
        "evidence_ms2_scan_number": pa.int64(),
        "evidence_ms2_native_id": pa.string(),
        "evidence_ms2_spectrum_index": pa.int64(),
        "fragment_mz": pa.float64(),
        "fragment_relative_intensity": pa.float64(),
        "isotope_state": pa.string(),
    }
    nested = {
        "ms2_scan_numbers": pa.list_(pa.int64()),
        "ms2_native_ids": pa.list_(pa.string()),
        "ms2_spectrum_indices": pa.list_(pa.int64()),
        "ion_relations": pa.list_(pa.struct(list(relation_fields.items()))),
    }
    schema = [pa.field("feature_id", pa.int64(), nullable=False)]
    for field in FEATURE_FIELDS:
        types = get_args(field.type) or (field.type,)
        dtype = nested.get(field.name)
        if dtype is None:
            dtype = (
                pa.float64()
                if float in types
                else pa.int64()
                if int in types
                else pa.bool_()
                if bool in types
                else pa.string()
            )
        schema.append(pa.field(field.name, dtype))
    return pa.schema(
        schema,
        metadata={
            b"masscubex.table_schema": str(TABLE_SCHEMA_VERSION).encode(),
            b"masscubex.retention_time_unit": b"second",
            b"masscubex.peak_area_unit": b"intensity * second",
            b"masscubex.feature_id": b"one-based row id within this input file",
            b"masscubex.output_precision": json.dumps(OUTPUT_PRECISION, sort_keys=True).encode(),
        },
    )


def _rounded_relations(relations, feature_id):
    """Copy sparse evidence dictionaries; never mutate scientific results."""
    result = []
    for relation in relations:
        if not RELATION_COLUMNS.issuperset(relation):
            raise ValueError(f"Unknown ion evidence fields in feature {feature_id}")
        rounded = dict(relation)
        for name, digits in RELATION_DECIMAL_PLACES.items():
            value = relation.get(name)
            if value is not None:
                if not math.isfinite(value):
                    raise ValueError(f"Non-finite ion evidence in feature {feature_id}")
                rounded[name] = round(float(value), digits)
        result.append(rounded)
    return result


def feature_rows(features):
    """Validate finite scalars/evidence before either writer publishes a table."""
    for index, feature in enumerate(features, 1):
        row = {name: getattr(feature, name) for name in FEATURE_COLUMNS[1:]}
        if any(row[name] is not None and not math.isfinite(row[name]) for name in FLOAT_COLUMNS):
            raise ValueError(f"Non-finite feature value (numerical overflow) in feature {index}")
        row["ion_relations"] = _rounded_relations(feature.ion_relations, index)
        yield {"feature_id": index, **row}


def _write_tsv(path, features):
    count = 0
    with open(path, "w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(FEATURE_COLUMNS)
        for row in feature_rows(features):
            writer.writerow(
                [
                    ""
                    if row[name] is None
                    else format(float(row[name]), f".{DECIMAL_PLACES[name]}f")
                    if name in FLOAT_COLUMNS
                    else json.dumps(row[name], separators=(",", ":"), allow_nan=False)
                    if name in NESTED_COLUMNS or name in BOOLEAN_COLUMNS
                    else row[name]
                    for name in FEATURE_COLUMNS
                ]
            )
            count += 1
    return count


def _parquet_batch(features, first_id, schema):
    """Build Arrow columns directly, avoiding an intermediate dictionary per row."""
    import pyarrow as pa
    import pyarrow.compute as pc

    columns = [pa.array(range(first_id, first_id + len(features)), type=pa.int64())]
    for field in schema:
        # Arrow exposes field.name through a property. Resolve it once per
        # column, rather than once for every row in the batch.
        name = field.name
        if name == "feature_id":
            continue
        if name == "ion_relations":
            values = [
                _rounded_relations(feature.ion_relations, index)
                for index, feature in enumerate(features, first_id)
            ]
        elif name in FLOAT_COLUMNS:
            digits = DECIMAL_PLACES[name]
            values = [
                None
                if (value := getattr(feature, name)) is None
                else round(float(value), digits)
                for feature in features
            ]
        else:
            values = [getattr(feature, name) for feature in features]
        column = pa.array(values, type=field.type)
        if name in FLOAT_COLUMNS and pc.all(pc.is_finite(column)).as_py() is False:
            raise ValueError(f"Non-finite feature value (numerical overflow) in {name}")
        columns.append(column)
    return pa.Table.from_arrays(columns, schema=schema)


def _write_parquet(path, features):
    import pyarrow as pa
    import pyarrow.parquet as pq

    schema = parquet_schema()
    categories = [field.name for field in schema if pa.types.is_string(field.type)]
    count = 0
    with pq.ParquetWriter(
        path, schema, compression="zstd", compression_level=3, use_dictionary=categories
    ) as writer:
        iterator = iter(features)
        while batch := list(islice(iterator, ROW_GROUP_SIZE)):
            writer.write_table(_parquet_batch(batch, count + 1, schema))
            count += len(batch)
    return count


def write_features(path, features, output_format):
    """Write a temporary table; the pipeline owns atomic publication and cleanup."""
    if output_format == "parquet":
        count = _write_parquet(path, features)
    elif output_format == "tsv":
        count = _write_tsv(path, features)
    else:
        raise ValueError("Output format must be parquet or tsv")
    # Windows' _commit (used by os.fsync) requires a writable descriptor.
    with open(path, "r+b" if os.name == "nt" else "rb") as stream:
        os.fsync(stream.fileno())
    return count


def iter_feature_rows(path):
    """Read typed rows from either format, streaming bounded Parquet row batches.

    Parquet relationship structs include null fields for evidence that is absent
    in that relationship kind; these are omitted here to match TSV dictionaries.
    """
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        import pyarrow.parquet as pq

        with pq.ParquetFile(path) as source:
            for batch in source.iter_batches(batch_size=ROW_GROUP_SIZE, use_threads=False):
                for row in batch.to_pylist():
                    for relation in row["ion_relations"]:
                        keep_null = (
                            {"parent_ion", "fragment_ion"}
                            if relation["relation"] in {"in_source_fragment", "neutral_loss"}
                            else set()
                        )
                        if relation.get("ms2_supported"):
                            keep_null.add("evidence_ms2_scan_number")
                        for name in list(relation):
                            if relation[name] is None and name not in keep_null:
                                del relation[name]
                    yield row
    else:
        integer_columns = {"feature_id"} | {
            f.name for f in ALL_FIELDS if int in (get_args(f.type) or (f.type,))
        }
        with path.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                for name, value in row.items():
                    row[name] = (
                        None
                        if value == ""
                        else float(value)
                        if name in LEGACY_FLOAT_COLUMNS
                        else int(value)
                        if name in integer_columns
                        else json.loads(value)
                        if name in NESTED_COLUMNS
                        or name in BOOLEAN_COLUMNS
                        or (name == "adduct" and (value.startswith('["') or value == "[]"))
                        else value
                    )
                yield row
