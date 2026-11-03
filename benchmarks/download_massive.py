"""Download a deterministic, stratified 50-file sample of MSV000086989.

Uses public original mzXML exports (40 biological, 6 blank, 4 QC) over FTPS.
Downloads remain outside version control. Run from the repository root:
  python benchmarks/download_massive.py --directory validation/MSV000086989
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
import argparse
import hashlib
import json
import re
import subprocess


PAGE = "https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=2f73277b4e034948acebfdf1edab17ed"


def spread(rows, count):
    rows = sorted(rows, key=lambda row: row["name"])
    if len(rows) < count:
        raise ValueError(f"Expected at least {count} files, found {len(rows)}")
    return [rows[round(i * (len(rows) - 1) / (count - 1))] for i in range(count)]


def download(row, directory):
    path = directory / row["name"]
    remote = f"/v03/MSV000086989/{row['relative_path']}/{row['name']}"
    url = "ftp://massive-ftp.ucsd.edu" + remote
    if not path.exists() or path.stat().st_size != row["size"]:
        temporary = Path(str(path) + ".download")
        subprocess.run(
            [
                "curl",
                "--ssl-reqd",
                "--fail",
                "--silent",
                "--show-error",
                "--disable-epsv",
                "--connect-timeout",
                "30",
                "--max-time",
                "600",
                "--retry",
                "2",
                "--retry-delay",
                "10",
                url,
                "-o",
                str(temporary),
            ],
            check=True,
        )
        if temporary.stat().st_size != row["size"]:
            raise ValueError(f"Unexpected download size: {path}")
        temporary.replace(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    result = {
        "name": row["name"],
        "bytes": row["size"],
        "sha256": digest,
        "url": url,
        "group": row["group"],
    }
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--listing", type=Path, help="Previously saved dataset row_data JSON")
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    if args.listing:
        rows = json.loads(args.listing.read_text())
    else:
        with urlopen(PAGE, timeout=60) as response:
            text = response.read().decode()
        match = re.search(r"var dataset_files = (.*);", text)
        rows = json.loads(match.group(1))["row_data"]
    rows = [row for row in rows if row["collection"] == "peak" and row["name"].endswith(".mzXML")]
    selected = []
    for prefix, count, group in [
        ("NIST_POS_Samp", 40, "biological"),
        ("Blank", 6, "blank"),
        ("QC_Mix", 4, "QC"),
    ]:
        selected.extend(
            {**row, "group": group}
            for row in spread([r for r in rows if r["name"].startswith(prefix)], count)
        )
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda row: download(row, args.directory), selected))
    (args.directory / "download_manifest.json").write_text(json.dumps(results, indent=2) + "\n")
    (args.directory / "files.txt").write_text("".join(row["name"] + "\n" for row in results))


if __name__ == "__main__":
    main()
