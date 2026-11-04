"""Lightweight input-format recognition for discovery and numerical readers."""


def file_format(path):
    name = str(path).lower()
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in ("mzml", "mzxml"):
        if name.endswith("." + suffix):
            return suffix
    raise ValueError(f"Expected .mzML, .mzXML, .mzML.gz, or .mzXML.gz: {path}")
