"""Lightweight input-format recognition for discovery and numerical readers."""

from contextlib import contextmanager
import gzip


@contextmanager
def open_xml(path):
    """Open plain/gzip XML by its bytes, including misnamed repository files."""
    with open(path, "rb") as stream:
        compressed = stream.read(2) == b"\x1f\x8b"
        stream.seek(0)
        if compressed:
            with gzip.GzipFile(fileobj=stream, mode="rb") as decoded:
                yield decoded
        else:
            yield stream


def file_format(path):
    name = str(path).lower()
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in ("mzml", "mzxml"):
        if name.endswith("." + suffix):
            return suffix
    raise ValueError(f"Expected .mzML, .mzXML, .mzML.gz, or .mzXML.gz: {path}")
