"""Initialize a file worker before allocating a large feature table.

Keep import/schema/provenance allocations out of the first large extraction.
The discarded one-row buffer only primes the existing writer; it is never
published, associated with input data or used to change scientific settings.
"""


def initialize(output_format="parquet"):
    from .pipeline import implementation_id
    from .features import warmup

    warmup()
    implementation_id()
    if output_format != "parquet":
        return
    from .features import Feature
    from .output import _write_parquet
    import pyarrow as pa

    feature = Feature(
        mz=100.0, rt_apex_seconds=1.0, rt_start_seconds=0.0, rt_end_seconds=2.0,
        peak_height=1.0, peak_area=1.0, scan_count=1, apex_noise_threshold=0.0,
        apex_threshold_ratio=None, polarity="positive",
    )
    sink = pa.BufferOutputStream()
    try:
        _write_parquet(sink, [feature])
    finally:
        sink.close()
