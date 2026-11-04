"""Export actual upstream trace samples and segmentation gates for EIC panels E/F.

PYTHONPATH=validation/output-v0.9/upstream/src:validation/output-v0.9/upstream-deps \
  .venv/bin/python benchmarks/export_upstream_eics.py
"""

import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from audit_upstream_features import mzml_records


def main():
    from masscube.params import Params
    from masscube.raw_data_utils.core import MSData
    from masscube.feature_evaluation import calculate_noise_score

    output = []
    for name, target in [
        ("Blank_00.mzXML", 171.1497344970703),
        ("QC_Mix_01.mzXML", 154.9901580810547),
    ]:
        source = Path("validation/MSV000086989") / name
        data = MSData()
        data.params = Params()
        data.params.set_default("orbitrap")
        data.update_metadata(
            {"file_path": str(source), "ms_type": "orbitrap", "ion_mode": "positive"}
        )
        data.extract_scans_mzml(mzml_records(source))
        data.detect_features()
        before = []
        for feature in data.features:
            if (
                abs(feature.mz - target) <= 0.01
                and feature.rt_seq[0] * 60 <= 100
                and feature.rt_seq[-1] * 60 >= 200
            ):
                score = float(calculate_noise_score(feature.signals[:, 1]))
                before.append(
                    dict(
                        mz=float(feature.mz),
                        noise_score=score,
                        skips_segmentation_because_noise_above_3=score > 3,
                    )
                )
        data.segment_features(iteration=2)
        data.finalize_features()
        archived = pq.read_table(
            Path("validation/output-v0.9/upstream-orbitrap") / (name + ".parquet")
        )
        assert len(data.features) == archived.num_rows
        np.testing.assert_array_equal(
            [float(f.mz) for f in data.features], archived["mz"].to_numpy()
        )
        np.testing.assert_array_equal(
            [float(f.rt) * 60 for f in data.features], archived["rt_apex_seconds"].to_numpy()
        )
        selected = []
        for i, feature in enumerate(data.features, 1):
            if (
                abs(feature.mz - target) <= 0.01
                and feature.rt_seq[0] * 60 <= 100
                and feature.rt_seq[-1] * 60 >= 200
            ):
                selected.append(
                    dict(
                        original_row=i,
                        mz=float(feature.mz),
                        time_seconds=(np.asarray(feature.rt_seq, dtype=float) * 60).tolist(),
                        intensity=feature.signals[:, 1].astype(float).tolist(),
                    )
                )
        result = dict(
            file=name,
            target_mz=target,
            original_features=len(data.features),
            pre_segmentation=before,
            final_traces=selected,
        )
        output.append(result)
        print(name, before, "exported traces:", len(selected), flush=True)
    Path("benchmarks/figures/eic-review-v0.9/upstream_traces.json").write_text(
        json.dumps(output, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
