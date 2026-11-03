"""Exercise optional grouping on a sample, blank and QC from the 50-file run.

Run validate_repository.py first. Verify that grouping preserves MS1 values,
MS2 assignments and ranking, and that every relationship is symmetric and
every claimed MS2 fragment is present in the independently decoded source.
Candidate counts are not validated compound identifications.
"""

from pathlib import Path
import argparse
import csv
import json
import time

import numpy as np

from masscubex import Config, __version__, process_file
from masscubex.features import warmup


def read_table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def verify_links(path, rows, config):
    from pyteomics import mzxml

    needed = {
        link["evidence_ms2_native_id"]
        for row in rows
        for link in json.loads(row["ion_relations"])
        if link.get("ms2_supported")
    }
    fragments = {}
    with mzxml.MzXML(str(path), use_index=False) as reader:
        for scan in reader:
            if str(scan["num"]) in needed:
                fragments[str(scan["num"])] = (scan["m/z array"], scan["intensity array"])
    assert set(fragments) == needed
    count = support = 0
    links = [json.loads(row["ion_relations"]) for row in rows]
    for index, row in enumerate(rows):
        for link in links[index]:
            count += 1
            partner = link["partner_feature_id"] - 1
            assert partner != index and 0 <= partner < len(rows)
            assert link["status"] == "candidate"
            assert row["polarity"] == rows[partner]["polarity"]
            assert link["shared_scans"] >= config.group_min_scans
            assert link["correlation"] >= config.group_min_correlation
            assert (
                abs(float(row["rt_apex_seconds"]) - float(rows[partner]["rt_apex_seconds"]))
                <= config.group_rt_tolerance + 1e-8
            )
            assert any(
                other["partner_feature_id"] == index + 1 and other["relation"] == link["relation"]
                for other in links[partner]
            )
            if link.get("ms2_supported"):
                support += 1
                parent = rows[link["parent_feature_id"] - 1]
                child = rows[link["fragment_feature_id"] - 1]
                assert link["evidence_ms2_native_id"] in json.loads(parent["ms2_native_ids"])
                assert link["evidence_ms2_scan_number"] == int(link["evidence_ms2_native_id"])
                mz, intensity = fragments[link["evidence_ms2_native_id"]]
                matches = np.flatnonzero(mz == link["fragment_mz"])
                assert len(matches)
                relative = float(intensity[matches[0]] / intensity.max())
                np.testing.assert_allclose(relative, link["fragment_relative_intensity"], rtol=1e-6)
                assert relative >= config.isf_min_relative_intensity - 1e-8
                assert (
                    abs(float(child["mz"]) - link["fragment_mz"])
                    <= max(
                        config.group_mz_tolerance_da,
                        link["fragment_mz"] * config.group_mz_tolerance_ppm * 1e-6,
                    )
                    + 1e-8
                )
    return {
        "symmetric_link_rows": count,
        "ms2_supported_link_rows": support,
        "evidence_ms2_spectra": len(needed),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument(
        "--output", type=Path, default=Path(f"benchmarks/results/ion-grouping-v{__version__}.json")
    )
    args = parser.parse_args()
    manifest = json.loads((args.directory / "download_manifest.json").read_text())
    selected = [
        next(row for row in manifest if row["group"] == group)
        for group in ("biological", "blank", "QC")
    ]
    config = Config(group_ions=True)
    warmup()
    results = []
    for row in selected:
        source = args.directory / row["name"]
        output = args.directory / f"v{__version__}" / "grouped" / (source.name + ".tsv")
        start = time.perf_counter()
        result = process_file(source, output, config=config, overwrite=True)
        elapsed = time.perf_counter() - start
        baseline = read_table(output.parent.parent / "features" / output.name)
        grouped = read_table(output)
        assert len(grouped) == len(baseline)
        for before, after in zip(baseline, grouped):
            assert {k: v for k, v in before.items() if k != "ion_relations"} == {
                k: v for k, v in after.items() if k != "ion_relations"
            }
        checked = verify_links(source, grouped, config)
        assert process_file(source, output, config=config, resume=True).status == "skipped"
        metadata = json.loads(Path(str(output) + ".json").read_text())
        results.append(
            {
                "name": source.name,
                "group": row["group"],
                "seconds": elapsed,
                "features": result.features,
                "ms1_scans": result.ms1_scans,
                **checked,
                **{
                    key: metadata[key]
                    for key in (
                        "adduct_relations",
                        "neutral_loss_relations",
                        "isf_relations",
                        "tsv_sha256",
                        "implementation_sha256",
                    )
                },
                "unchanged_measurements_and_ms2": True,
                "resume_verified": True,
            }
        )
        print(json.dumps(results[-1]), flush=True)
    report = {
        "version": __version__,
        "config": config.to_dict(),
        "results": results,
        "limitation": "Checks verify candidate evidence and bookkeeping, not chemical identity or false-positive rates.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
