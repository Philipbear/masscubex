"""Render explicitly selected EIC audit cases from archived, unchanged outputs."""

import json
import os
from pathlib import Path

import numpy as np

from audit_eics import extract_eics, raw_ms1, tables

os.environ.setdefault("MPLCONFIGDIR", "/tmp/masscubex-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import EngFormatter


ROOT = Path("validation")
OUT = Path("benchmarks/figures/eic-review-v0.9")
CASES = [
    ("Blank_00.mzXML", 173.07944454017115, 8, 48, "A  Blank: small peaks amid background"),
    ("QC_Mix_01.mzXML", 248.15050130727388, 195, 244, "B  QC: a bounded weak peak"),
    ("Blank_Post_10-02.mzXML", 199.18134, 0, 600, "C  Blank: persistent low-level signal"),
    ("QC_Mix_01.mzXML", 199.18134, 0, 600, "D  QC: persistent low-level signal"),
    ("Blank_00.mzXML", 171.1497344970703, 100, 200, "E  Blank: splitting persistent background"),
    ("QC_Mix_01.mzXML", 154.9901580810547, 100, 200, "F  QC: splitting persistent background"),
]


def matches(table, mass, radius, start, end):
    return table[
        (abs(table.mz - mass) <= radius)
        & (table.rt_end_seconds >= start)
        & (table.rt_start_seconds <= end)
    ]


def plot():
    OUT.mkdir(parents=True, exist_ok=True)
    original_traces = json.loads((OUT / "upstream_traces.json").read_text())
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    report = []
    for ax, (name, mass, start, end, title) in zip(axes.flat, CASES):
        current, original, floor, config = tables(ROOT, name)
        t, masses, y, offsets = raw_ms1(ROOT / "MSV000086989" / name, ROOT / "eic-review/raw")
        radius = max(config["mz_tolerance_da"], mass * config["mz_tolerance_ppm"] * 1e-6)
        narrow, broad = extract_eics(
            masses, y, offsets, np.array([mass, mass]), np.array([radius, 0.01])
        )
        view = (t >= start) & (t <= end)
        ax.plot(t[view], broad[view], color="#b8c2cd", lw=1.2, alpha=0.8)
        ax.plot(t[view], narrow[view], color="#182e40", lw=1, marker=".", markersize=2)
        for entry in original_traces:
            if entry["file"] == name and abs(entry["target_mz"] - mass) < radius:
                for trace in entry["final_traces"]:
                    original_t = np.array(trace["time_seconds"])
                    original_y = np.array(trace["intensity"])
                    visible = (original_t >= start) & (original_t <= end)
                    ax.plot(
                        original_t[visible],
                        original_y[visible],
                        color="#197b9d",
                        lw=0.9,
                        ls="--",
                        alpha=0.85,
                    )
        if narrow[view].max() <= 50000:
            ax.axhline(50000, color="#8b7145", lw=1, ls="--")
        counts = {}
        for label, table, tolerance, ypos, color in [
            ("MassCubeX", current, radius, -0.13, "#d35d32"),
            ("X, 50k floor", floor, radius, -0.22, "#9766b1"),
            ("Original", original, 0.01, -0.31, "#197b9d"),
        ]:
            selected = matches(table, mass, tolerance, start, end)
            apex_count = int(selected.rt_apex_seconds.between(start, end).sum())
            counts[label] = dict(
                overlapping_rows=len(selected),
                apices_in_view=apex_count,
                apices_whole_run=int((abs(table.mz - mass) <= tolerance).sum()),
            )
            for row in selected.itertuples():
                ax.plot(
                    [max(start, row.rt_start_seconds), min(end, row.rt_end_seconds)],
                    [ypos, ypos],
                    color=color,
                    lw=2,
                    alpha=0.65,
                    transform=ax.get_xaxis_transform(),
                    clip_on=False,
                )
                if start <= row.rt_apex_seconds <= end:
                    ax.plot(
                        row.rt_apex_seconds,
                        ypos,
                        "|",
                        color=color,
                        markersize=8,
                        transform=ax.get_xaxis_transform(),
                        clip_on=False,
                    )
            ax.text(
                0,
                ypos - 0.035,
                f"{label}: {len(selected)} intervals ({apex_count} apices)",
                transform=ax.transAxes,
                fontsize=8,
                color=color,
                va="top",
            )
        ax.set_title(f"{title}\n{name} · m/z {mass:.5f}", loc="left", fontsize=10.5)
        ax.set_xlim(start, end)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_formatter(EngFormatter(sep=""))
        ax.set_ylabel("Raw intensity")
        ax.set_xlabel("Retention time (s)", labelpad=43)
        ax.grid(axis="y", alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
        report.append(
            dict(
                file=name,
                mz=mass,
                view_seconds=[start, end],
                raw_eic_radius_da=radius,
                counts=counts,
                raw_max_in_view=float(narrow[view].max()),
                raw_median_in_view=float(np.median(narrow[view])),
                nonzero_fraction_in_view=float(np.mean(narrow[view] > 0)),
            )
        )
    fig.suptitle(
        "Extra features before algorithm changes: raw EICs and reported peak intervals", fontsize=15
    )
    fig.legend(
        handles=[
            Line2D([], [], color="#182e40", label="Raw EIC: learned MS1 tolerance"),
            Line2D([], [], color="#b8c2cd", label="Raw EIC: ±0.01 Da for context"),
            Line2D([], [], color="#8b7145", ls="--", label="Original 50,000 intensity floor"),
            Line2D([], [], color="#197b9d", ls="--", label="Actual original trace (E/F)"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.subplots_adjust(top=0.89, bottom=0.08, hspace=0.95, wspace=0.2)
    fig.savefig(OUT / "eic_comparison.png", dpi=180)
    fig.savefig(OUT / "eic_comparison.pdf")
    (OUT / "eic_comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    # Two-panel views remain legible when embedded in a chat.
    for row, name in enumerate(("weak_peaks", "persistent_signal", "splitting")):
        fig.canvas.draw()
        boxes = [
            ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
            for ax in axes[row]
        ]
        from matplotlib.transforms import Bbox

        fig.savefig(
            OUT / (name + ".png"), dpi=180, bbox_inches=Bbox.union(boxes).expanded(1.02, 1.08)
        )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    plot()
