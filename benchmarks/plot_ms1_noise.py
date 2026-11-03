"""Plot real MS1 centroids before/after the implemented adaptive filter.

Requires matplotlib in addition to the runtime dependencies. By default use
the first sample file, the scan nearest the midpoint of acquisition time, and
a fixed m/z 100–250 window. Thresholds are computed from the complete scan.
This diagnostic does not change processing defaults or label peaks as truth.
"""

from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import tempfile

import numpy as np

from masscubex import Config, __version__
from masscubex.noise import estimate_noise
from masscubex.pipeline import implementation_id
from masscubex.readers import iter_ms1


def select_scan(path, target_time):
    times = np.array([scan.rt for scan in iter_ms1(path)])
    if not len(times):
        raise ValueError("No MS1 scans in input")
    requested = float((times[0] + times[-1]) / 2) if target_time is None else target_time
    ordinal = int(np.argmin(abs(times - requested)))
    scans = iter_ms1(path)
    try:
        for index, scan in enumerate(scans):
            if index == ordinal:
                return scan, ordinal, requested, (float(times[0]), float(times[-1]))
    finally:
        scans.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=Path("validation/MSV000086989/NIST_POS_Samp_01-01.mzXML"),
    )
    parser.add_argument("--rt-seconds", type=float)
    parser.add_argument("--mz-min", type=float, default=100)
    parser.add_argument("--mz-max", type=float, default=250)
    parser.add_argument("--noise-sigma", type=float, choices=(3.0,), default=3)
    parser.add_argument("--compare-sigma", type=float, default=2)
    parser.add_argument("--axis-multiple", type=float, default=3.5)
    parser.add_argument("--out-dir", type=Path, default=Path("benchmarks/figures/legacy_clipped"))
    args = parser.parse_args()
    if not 0 < args.compare_sigma < args.noise_sigma or args.axis_multiple <= 1:
        parser.error("Require 0 < compare-sigma < noise-sigma and axis-multiple > 1")
    if not args.mz_min < args.mz_max:
        parser.error("Require mz-min < mz-max")
    config = Config(
        noise_estimator="clipped",
        noise_sigma=args.noise_sigma,
        noise_window_da=25,
        noise_min_neighbors=8,
    )
    scan, ordinal, requested, time_range = select_scan(args.input, args.rt_seconds)
    if not scan.centroided:
        raise ValueError("This before/after centroid diagnostic expects centroided MS1 input")
    background, std, cutoff = estimate_noise(scan.mz, scan.intensity, config)
    # Historical plotting comparison only. Production Config fixes sigma at 3.
    # The background fit does not depend on the detection multiplier.
    relaxed = np.maximum(config.min_intensity, background + args.compare_sigma * std)
    kept = scan.intensity > cutoff
    kept_relaxed = scan.intensity > relaxed
    additional = kept_relaxed & ~kept
    window = (scan.mz >= args.mz_min) & (scan.mz <= args.mz_max)
    if not window.any() or np.median(background[window]) <= 0 or np.median(cutoff[window]) <= 0:
        raise ValueError("Choose a window with supported, nonzero background estimates")
    assert np.all(kept <= kept_relaxed)
    median_background = float(np.median(background[window]))
    median_cutoff = float(np.median(cutoff[window]))
    limits = [args.axis_multiple * median_cutoff, args.axis_multiple * median_background]

    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "masscubex-mpl"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#aeb8c2",
            "axes.labelcolor": "#243548",
            "text.color": "#243548",
            "xtick.color": "#526170",
            "ytick.color": "#526170",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
        }
    )
    colors = {
        "raw": "#657785",
        "kept": "#087f8c",
        "extra": "#d38912",
        "mean": "#8c63a5",
        "cutoff": "#d64c4c",
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)

    def number(value, position=None):
        return "0" if value == 0 else f"{value / 1000:,.0f}k"

    def sticks(ax, selection, color, limit, width=1.1):
        x, y = scan.mz[selection], scan.intensity[selection]
        ax.vlines(x, 0, y, color=color, linewidth=width, alpha=0.93, zorder=2)
        over = y > limit
        ax.scatter(
            x[over],
            np.full(over.sum(), limit * 0.975),
            marker="^",
            s=22,
            color=color,
            linewidths=0,
            zorder=5,
        )
        return int(over.sum())

    def panel(ax, row, selection, limit, lo, hi):
        active_cutoff = relaxed if row == 2 else cutoff
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, limit)
        ax.grid(axis="y", color="#e8edf1", linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.yaxis.set_major_formatter(FuncFormatter(number))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        n_clipped = 0
        if row == 0:
            n_clipped = sticks(ax, selection, colors["raw"], limit)
        else:
            n_clipped = sticks(ax, selection & kept, colors["kept"], limit, 1.5)
            if row == 2:
                n_clipped += sticks(ax, selection & additional, colors["extra"], limit, 1.7)
        ax.plot(
            scan.mz[selection], background[selection], color=colors["mean"], linewidth=1.3, zorder=3
        )
        ax.plot(
            scan.mz[selection],
            active_cutoff[selection],
            color=colors["cutoff"],
            linewidth=1.4,
            linestyle=(0, (4, 2.5)),
            zorder=4,
        )
        cutoff_hidden = int(np.count_nonzero(selection & (active_cutoff > limit)))
        label = f"{n_clipped} tall peaks clipped"
        if cutoff_hidden:
            label += f"  |  cutoff off-scale at {cutoff_hidden} points"
        ax.text(
            0.985,
            0.92,
            label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 2},
        )
        ax.set_ylabel("Intensity (a.u.)")

    n = int(window.sum())
    n_default, n_relaxed = (
        int(np.count_nonzero(window & kept)),
        int(np.count_nonzero(window & kept_relaxed)),
    )
    fig, axes = plt.subplots(3, 2, figsize=(14, 10.5))
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.17, top=0.80, hspace=0.39, wspace=0.22)
    fig.suptitle(
        "Previous clipped MS1 filter: before / after",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=21,
        weight="bold",
    )
    fig.text(
        0.075,
        0.936,
        f"{args.input.name}  •  scan {scan.native_id}  •  RT {scan.rt:.3f} s  •  m/z {args.mz_min:g}–{args.mz_max:g}",
        fontsize=12,
        color="#526170",
    )
    for x, heading, description in (
        (
            0.075,
            "Cutoff view",
            f"Median clipped cutoff = {median_cutoff / 1000:.1f}k; {100 / args.axis_multiple:.0f}% of axis height",
        ),
        (
            0.575,
            "Background zoom",
            f"Median estimated background = {median_background / 1000:.1f}k; {100 / args.axis_multiple:.0f}% of axis height",
        ),
    ):
        fig.text(x, 0.885, heading, fontsize=16, weight="bold")
        fig.text(x, 0.856, description, fontsize=10.5, color="#526170")
    titles = (
        f"Before filtering  ·  {n} centroids",
        f"Previous clipped {args.noise_sigma:g}σ  ·  {n_default}/{n} kept ({n_default / n:.1%})",
        f"Less strict {args.compare_sigma:g}σ  ·  {n_relaxed}/{n} kept (+{n_relaxed - n_default})",
    )
    for row in range(3):
        for col in range(2):
            ax = axes[row, col]
            panel(ax, row, window, limits[col], args.mz_min, args.mz_max)
            ax.set_title(titles[row], loc="left", pad=10, weight="bold")
            if row == 2:
                ax.set_xlabel("m/z")
    handles = [
        Line2D([], [], color=colors["raw"], label="Original centroids"),
        Line2D([], [], color=colors["kept"], label=f"Retained at {args.noise_sigma:g}σ"),
        Line2D([], [], color=colors["extra"], label=f"Additional at {args.compare_sigma:g}σ"),
        Line2D([], [], color=colors["mean"], label="Estimated background μ"),
        Line2D([], [], color=colors["cutoff"], linestyle="--", label="Active cutoff μ + kσ"),
    ]
    fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.065, 0.075),
        ncol=3,
        frameon=False,
        fontsize=10.5,
    )
    fig.text(
        0.075,
        0.055,
        "Same raw intensities in both columns; only the y-axis limits change. ▲ marks peaks extending above the axis.",
        fontsize=10,
        color="#526170",
    )
    fig.text(
        0.075,
        0.032,
        "Thresholds use the complete scan: 25 Da window, 2 clipping passes at 2σ, minimum 8 neighbors. Removed ≠ confirmed noise.",
        fontsize=10,
        color="#526170",
    )
    fig.savefig(args.out_dir / "ms1_noise_before_after.png", dpi=160)
    fig.savefig(args.out_dir / "ms1_noise_before_after.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(14, 9))
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.16, top=0.86, hspace=0.42)
    fig.suptitle(
        "Full-scan context · same MS1 acquisition",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=20,
        weight="bold",
    )
    fig.text(
        0.075,
        0.932,
        f"{args.input.name}  •  scan {scan.native_id}  •  RT {scan.rt:.3f} s  •  linear intensity scale",
        fontsize=12,
    )
    all_points = np.ones(len(scan.mz), dtype=bool)
    full_titles = (
        f"Before filtering · {len(scan.mz)} centroids",
        f"Previous clipped {args.noise_sigma:g}σ · {int(kept.sum())}/{len(scan.mz)} kept ({kept.mean():.1%})",
        f"Less strict {args.compare_sigma:g}σ · {int(kept_relaxed.sum())}/{len(scan.mz)} kept ({kept_relaxed.mean():.1%})",
    )
    for row, ax in enumerate(axes):
        panel(ax, row, all_points, limits[0], float(scan.mz[0]), float(scan.mz[-1]))
        ax.axvspan(args.mz_min, args.mz_max, color="#dceaf1", alpha=0.45, zorder=0)
        ax.set_title(full_titles[row], loc="left", weight="bold", pad=9)
    axes[-1].set_xlabel("m/z")
    fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.065, 0.055),
        ncol=3,
        frameon=False,
        fontsize=10.5,
    )
    fig.text(
        0.075,
        0.025,
        "Shaded band = detailed view. Sparse windows use the zero floor; a zero cutoff means insufficient local background support.",
        fontsize=10,
        color="#526170",
    )
    fig.savefig(args.out_dir / "ms1_noise_full_scan.png", dpi=160)
    plt.close(fig)

    with (args.out_dir / "ms1_noise_scan.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "mz",
                "original_intensity",
                "estimated_background",
                "background_sd",
                "default_cutoff",
                "comparison_cutoff",
                "kept_default",
                "kept_comparison",
            ]
        )
        writer.writerows(
            zip(
                scan.mz,
                scan.intensity,
                background,
                std,
                cutoff,
                relaxed,
                kept.astype(int),
                kept_relaxed.astype(int),
            )
        )
    counts = {}
    for name, mask in (("shown_window", window), ("whole_scan", all_points)):
        total_intensity = float(scan.intensity[mask].sum())
        counts[name] = {
            "input_centroids": int(mask.sum()),
            "kept_default": int(np.count_nonzero(mask & kept)),
            "kept_comparison": int(np.count_nonzero(mask & kept_relaxed)),
            "removed_default": int(np.count_nonzero(mask & ~kept)),
            "additional_at_comparison": int(np.count_nonzero(mask & additional)),
            "intensity_fraction_retained_default": float(
                scan.intensity[mask & kept].sum() / total_intensity
            ),
            "intensity_fraction_retained_comparison": float(
                scan.intensity[mask & kept_relaxed].sum() / total_intensity
            ),
            "zero_cutoff_points": int(np.count_nonzero(mask & (cutoff == 0))),
        }
    metadata = {
        "version": __version__,
        "source_file": args.input.name,
        "source_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "selection": {
            "scan_rule": "nearest_acquisition_time_midpoint"
            if args.rt_seconds is None
            else "nearest_requested_retention_time",
            "mass_window_rule": "fixed_requested_range",
            "default_window": [100, 250],
        },
        "requested_rt_seconds": requested,
        "acquisition_rt_range_seconds": time_range,
        "native_id": scan.native_id,
        "spectrum_index": scan.spectrum_index,
        "ms1_ordinal_zero_based": ordinal,
        "rt_seconds": scan.rt,
        "centroided": scan.centroided,
        "polarity": scan.polarity,
        "mz_window": [args.mz_min, args.mz_max],
        "config": config.to_dict(),
        "comparison_noise_sigma": args.compare_sigma,
        "implementation_sha256": implementation_id()[0],
        "median_background_in_window": median_background,
        "median_default_cutoff_in_window": median_cutoff,
        "axis_multiple": args.axis_multiple,
        "ymax_cutoff_view": limits[0],
        "ymax_background_view": limits[1],
        "counts": counts,
        "interpretation": "These are centroids removed by a rule, not known noise labels. The figure is not ground-truth validation or a change of defaults.",
    }
    (args.out_dir / "ms1_noise_example.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(
        json.dumps(
            {
                "output_directory": str(args.out_dir.resolve()),
                "native_id": scan.native_id,
                "rt_seconds": scan.rt,
                "counts": counts,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
