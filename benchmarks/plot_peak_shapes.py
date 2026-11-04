"""Before/after intervals on unchanged raw EICs, including original MassCube."""

import json
import os
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from audit_eics import extract_eics, raw_ms1
from plot_eic_review import CASES, matches

os.environ.setdefault("MPLCONFIGDIR", "/tmp/masscubex-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import EngFormatter
from matplotlib.transforms import Bbox

ROOT = Path("validation")
OUT = Path("benchmarks/figures/peak-shapes-v0.11")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    original_traces = json.loads(
        Path("benchmarks/figures/eic-review-v0.9/upstream_traces.json").read_text()
    )
    fig, axes = plt.subplots(3, 2, figsize=(14, 13))
    report = []
    for ax, (name, mass, start, end, title) in zip(axes.flat, CASES):
        before_path = ROOT / "eic-review/v0.10/parquet" / (name + ".parquet")
        after_path = ROOT / "peak-shape-v0.11/parquet" / (name + ".parquet")
        before = pq.read_table(before_path).to_pandas()
        after = pq.read_table(after_path).to_pandas()
        before = before[before.feature_origin == "ms1_detected"]
        after = after[after.feature_origin == "ms1_detected"]
        original = pq.read_table(
            ROOT / "output-v0.9/upstream-orbitrap" / (name + ".parquet")
        ).to_pandas()
        config = json.loads(Path(str(after_path) + ".json").read_text())["parameter_tuning"][
            "effective_configs"
        ]["positive"]
        radius = max(config["mz_tolerance_da"], mass * config["mz_tolerance_ppm"] * 1e-6)
        t, masses, intensity, offsets = raw_ms1(
            ROOT / "MSV000086989" / name, ROOT / "eic-review/raw"
        )
        raw = extract_eics(masses, intensity, offsets, np.array([mass]), np.array([radius]))[0]
        view = (t >= start) & (t <= end)
        ax.plot(t[view], raw[view], color="#243b4b", lw=1, marker=".", markersize=2)
        for entry in original_traces:
            if entry["file"] == name and abs(entry["target_mz"] - mass) < radius:
                for trace in entry["final_traces"]:
                    rt, y = np.array(trace["time_seconds"]), np.array(trace["intensity"])
                    keep = (rt >= start) & (rt <= end)
                    ax.plot(rt[keep], y[keep], color="#3a83ac", lw=0.9, ls="--", alpha=0.65)
        counts = {}
        for label, table, tolerance, position, color in [
            ("Before (v0.10)", before, radius, -0.13, "#d16d45"),
            ("After (v0.11)", after, radius, -0.22, "#168579"),
            ("Original MassCube", original, 0.01, -0.31, "#3a83ac"),
        ]:
            selected = matches(table, mass, tolerance, start, end)
            apex_count = int(selected.rt_apex_seconds.between(start, end).sum())
            counts[label] = dict(
                overlapping_intervals=len(selected),
                apices_in_view=apex_count,
                intervals_whole_run=int((abs(table.mz - mass) <= tolerance).sum()),
            )
            for row in selected.itertuples():
                ax.plot(
                    [max(start, row.rt_start_seconds), min(end, row.rt_end_seconds)],
                    [position, position],
                    color=color,
                    lw=2.4,
                    alpha=0.8,
                    transform=ax.get_xaxis_transform(),
                    clip_on=False,
                )
                if start <= row.rt_apex_seconds <= end:
                    ax.plot(
                        row.rt_apex_seconds,
                        position,
                        "|",
                        color=color,
                        markersize=9,
                        transform=ax.get_xaxis_transform(),
                        clip_on=False,
                    )
            interval_label = "interval" if len(selected) == 1 else "intervals"
            apex_label = "apex" if apex_count == 1 else "apices"
            ax.text(
                0,
                position - 0.035,
                f"{label}: {len(selected)} {interval_label} ({apex_count} {apex_label})",
                transform=ax.transAxes,
                fontsize=8.5,
                color=color,
                va="top",
            )
        ax.set_title(f"{title}\n{name} · m/z {mass:.5f}", loc="left", fontsize=10.5)
        ax.set_xlim(start, end)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("Raw MS1 intensity")
        ax.set_xlabel("Retention time (s)", labelpad=45)
        ax.yaxis.set_major_formatter(EngFormatter(sep=""))
        ax.grid(axis="y", alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
        entry = dict(
            file=name,
            mz=mass,
            view_seconds=[start, end],
            raw_eic_radius_da=radius,
            counts=counts,
            after_features=[],
        )
        for row in matches(after, mass, radius, start, end).itertuples():
            entry["after_features"].append(
                dict(
                    feature_id=int(row.feature_id),
                    mz=float(row.mz),
                    start_seconds=float(row.rt_start_seconds),
                    apex_seconds=float(row.rt_apex_seconds),
                    end_seconds=float(row.rt_end_seconds),
                    height=float(row.peak_height),
                    gaussian_similarity=float(row.gaussian_similarity),
                    segmentation_status=row.segmentation_status,
                )
            )
        # The bounded QC example is an independent preservation control.
        if name == "QC_Mix_01.mzXML" and abs(mass - 248.15050) < 1e-4:
            columns = [
                "mz",
                "rt_apex_seconds",
                "rt_start_seconds",
                "rt_end_seconds",
                "peak_height",
                "peak_area",
                "scan_count",
            ]
            old = matches(before, mass, radius, start, end)[columns].to_dict("records")
            new = matches(after, mass, radius, start, end)[columns].to_dict("records")
            assert len(old) == len(new) == 1 and old == new
            entry["weak_qc_measurements_unchanged"] = True
        report.append(entry)
    fig.suptitle("Shape guard: fewer background intervals, unchanged raw MS1 signals", fontsize=15)
    fig.legend(
        handles=[
            Line2D([], [], color="#243b4b", label="Raw EIC at the learned MS1 tolerance"),
            Line2D([], [], color="#3a83ac", ls="--", label="Actual original MassCube trace (E/F)"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.953),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.text(
        0.5,
        0.016,
        "Bars = reported intervals; ticks = raw apices. Remaining intervals can reflect acquisition/mass-trace gaps.\nCounts alone do not establish contamination or chemical peak identity.",
        ha="center",
        fontsize=9,
        color="#455562",
    )
    fig.subplots_adjust(top=0.86, bottom=0.14, hspace=0.94, wspace=0.20)
    fig.savefig(OUT / "comparison.png", dpi=180)
    # Crops contain only their axes: figure-wide captions/legends otherwise
    # intrude into the first/last panel row after the tight bounding box shifts.
    for artist in [*fig.texts, *fig.legends]:
        artist.set_visible(False)
    for index, name in enumerate(["weak_peaks", "persistent_signal", "splitting"]):
        fig.canvas.draw()
        boxes = [
            ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
            for ax in axes[index]
        ]
        fig.savefig(
            OUT / (name + ".png"), dpi=180, bbox_inches=Bbox.union(boxes).expanded(1.02, 1.08)
        )
    (OUT / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps([dict(file=r["file"], mz=r["mz"], counts=r["counts"]) for r in report], indent=2)
    )


if __name__ == "__main__":
    main()
