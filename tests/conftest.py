import base64
import gzip
import zlib

import numpy as np
import pytest


def encode(values, precision=64, endian="<", compress=True):
    raw = np.asarray(values, dtype=endian + ("f8" if precision == 64 else "f4")).tobytes()
    return base64.b64encode(zlib.compress(raw) if compress else raw).decode()


def mzml_spectrum(
    index,
    rt,
    mz,
    intensity,
    *,
    level=1,
    mode="centroid",
    unit="UO:0000010",
    polarity="positive",
    precision=64,
    compress=True,
    refs=False,
    precursor="",
):
    def cv(accession, value="", **attrs):
        extra = " ".join(f'{k}="{v}"' for k, v in attrs.items())
        return f'<cvParam accession="{accession}" value="{value}" {extra}/>'

    params = cv("MS:1000511", level)
    if mode:
        params += cv("MS:1000127" if mode == "centroid" else "MS:1000128")
    if polarity:
        params += cv(
            "MS:1000130"
            if polarity == "positive"
            else "MS:1000139"
            if polarity == "bogus"
            else "MS:1000129"
        )
    if refs:
        params = '<referenceableParamGroupRef ref="spectrum_params"/>'
    arrays = ""
    for key, values in [("MS:1000514", mz), ("MS:1000515", intensity)]:
        terms = cv(key) + cv("MS:1000523" if precision == 64 else "MS:1000521")
        terms += cv("MS:1000574" if compress else "MS:1000576")
        arrays += f"<binaryDataArray>{terms}<binary>{encode(values, precision, compress=compress)}</binary></binaryDataArray>"
    time = cv("MS:1000016", rt, unitAccession=unit)
    return (
        f'<spectrum id="scan={index}" defaultArrayLength="{len(mz)}">{params}'
        f"<scanList><scan>{time}</scan></scanList>"
        f"{precursor}"
        f"<binaryDataArrayList>{arrays}</binaryDataArrayList></spectrum>"
    )


def write_mzml(path, spectra, *, groups="", indexed=False):
    text = (
        '<mzML xmlns="http://psi.hupo.org/ms/mzml">'
        + groups
        + "<run><spectrumList>"
        + "".join(spectra)
        + "</spectrumList></run></mzML>"
    )
    if indexed:
        text = '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">' + text + "</indexedmzML>"
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wb") as stream:
        stream.write(text.encode())
    return path


def mzxml_scan(
    index,
    rt,
    mz,
    intensity,
    *,
    level=1,
    centroid="1",
    polarity="+",
    precision=64,
    endian=">",
    compress=True,
    nested="",
    precursor="",
):
    pairs = np.column_stack((mz, intensity)).reshape(-1)
    encoded = encode(pairs, precision, endian, compress)
    mode = f'centroided="{centroid}"' if centroid is not None else ""
    order = "network" if endian == ">" else "little"
    compression = "zlib" if compress else "none"
    return (
        f'<scan num="{index}" msLevel="{level}" peaksCount="{len(mz)}" {mode} '
        f'polarity="{polarity}" retentionTime="PT{rt}S">'
        f"{precursor}"
        f'<peaks precision="{precision}" byteOrder="{order}" compressionType="{compression}" '
        f'pairOrder="m/z-int">{encoded}</peaks>{nested}</scan>'
    )


def write_mzxml(path, scans, *, processing=""):
    text = f'<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2"><msRun>{processing}{"".join(scans)}</msRun></mzXML>'
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wb") as stream:
        stream.write(text.encode())
    return path


@pytest.fixture
def sample_file(tmp_path):
    y = [10, 20, 40, 80, 100, 80, 40, 20, 10]
    return write_mzml(
        tmp_path / "sample.mzML",
        [mzml_spectrum(i, i * 2, [100.123456789], [height]) for i, height in enumerate(y)],
    )
