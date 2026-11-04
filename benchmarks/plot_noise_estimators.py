"""Compare historical MS1 settings and the current low-tail default.

Use the same fixed real scan and m/z range as plot_ms1_noise.py. All thresholds
come from the full scan; plotting never changes data or processing settings.
Requires matplotlib (a diagnostic dependency, not a runtime dependency).
"""

from dataclasses import replace
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
from plot_ms1_noise import select_scan


def main():
    defaults = Config()
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
    parser.add_argument("--noise-fraction", type=float, default=defaults.noise_fraction)
    parser.add_argument("--noise-min-neighbors", type=int, default=defaults.noise_min_neighbors)
    parser.add_argument("--noise-window-da", type=float, default=defaults.noise_window_da)
    parser.add_argument(
        "--noise-layout", choices=("bins", "sliding"), default=defaults.noise_layout
    )
    parser.add_argument("--noise-sigma", type=float, choices=(3.0,), default=defaults.noise_sigma)
    parser.add_argument(
        "--comparison",
        choices=("five-percent", "previous-defaults", "clipped", "v0.9"),
        default="five-percent",
    )
    parser.add_argument("--axis-multiple", type=float, default=3.5)
    parser.add_argument("--out-dir", type=Path, default=Path(f"benchmarks/figures/v{__version__}"))
    args = parser.parse_args()
    if not args.mz_min < args.mz_max or args.axis_multiple <= 1:
        parser.error("Require mz-min < mz-max and axis-multiple > 1")
    current = Config(
        noise_estimator="low_tail",
        noise_fraction=args.noise_fraction,
        noise_min_neighbors=args.noise_min_neighbors,
        noise_window_da=args.noise_window_da,
        noise_layout=args.noise_layout,
        noise_sigma=args.noise_sigma,
    )
    if args.comparison == "five-percent":
        previous = replace(current, noise_fraction=0.05)
        previous_label = f"Lowest 5%, min {previous.noise_min_neighbors}"
    elif args.comparison == "v0.9":
        previous = replace(current, noise_layout="sliding", noise_window_da=50)
        previous_label = "v0.9: lowest 10%, min 5, sliding"
    else:
        previous = Config(
            noise_layout="sliding",
            noise_estimator="clipped" if args.comparison == "clipped" else "low_tail",
            noise_window_da=25,
            noise_fraction=0.05,
            noise_min_neighbors=8,
            noise_sigma=3,
        )
        previous_label = (
            "Previous clipped" if args.comparison == "clipped" else "Previous 5%, min 8"
        )
    scan, ordinal, requested, time_range = select_scan(args.input, args.rt_seconds)
    if not scan.centroided:
        raise ValueError("This diagnostic expects a centroid MS1 scan")
    window = (scan.mz >= args.mz_min) & (scan.mz <= args.mz_max)
    all_points = np.ones(len(scan.mz), dtype=bool)
    new_mean, new_sd, new_cutoff = estimate_noise(scan.mz, scan.intensity, current)
    old_mean, old_sd, old_cutoff = estimate_noise(scan.mz, scan.intensity, previous)
    new_keep, old_keep = scan.intensity > new_cutoff, scan.intensity > old_cutoff
    extra = new_keep & ~old_keep
    positive_mz = scan.mz[scan.intensity > 0]
    if current.noise_layout == "bins":
        edges = np.floor(scan.mz / current.noise_window_da) * current.noise_window_da
        neighbors = np.searchsorted(positive_mz, edges + current.noise_window_da) - np.searchsorted(
            positive_mz, edges
        )
    else:
        half = current.noise_window_da / 2
        neighbors = (
            np.searchsorted(positive_mz, scan.mz + half, side="right")
            - np.searchsorted(positive_mz, scan.mz - half, side="left")
            - (scan.intensity > 0).astype(int)
        )
    selected = np.where(
        neighbors >= current.noise_min_neighbors,
        np.maximum(current.noise_min_neighbors, np.ceil(neighbors * current.noise_fraction)),
        0,
    ).astype(int)
    supported = selected > 0
    if not window.any() or np.median(new_mean[window]) <= 0:
        raise ValueError("Choose a window with supported, nonzero background")
    limits = [
        args.axis_multiple * float(np.median(old_cutoff[window])),
        args.axis_multiple * float(np.median(new_cutoff[window])),
    ]
    counts = {}
    for name, mask in (("shown_window", window), ("whole_scan", all_points)):
        total_intensity = scan.intensity[mask].sum()
        valid = mask & supported
        counts[name] = {
            "input_centroids": int(mask.sum()),
            "kept_previous": int((mask & old_keep).sum()),
            "kept_low_tail": int((mask & new_keep).sum()),
            "additional_low_tail": int((mask & extra).sum()),
            "lost_previous": int((mask & old_keep & ~new_keep).sum()),
            "intensity_fraction_previous": float(
                scan.intensity[mask & old_keep].sum() / total_intensity
            ),
            "intensity_fraction_low_tail": float(
                scan.intensity[mask & new_keep].sum() / total_intensity
            ),
            "unsupported_low_tail": int((mask & ~supported).sum()),
            "minimum_expands_fraction": int(
                (valid & (np.ceil(neighbors * current.noise_fraction) < selected)).sum()
            ),
            "selected_count_range_supported": [
                int(selected[valid].min()),
                int(selected[valid].max()),
            ],
            "effective_fraction_range_supported": [
                float((selected[valid] / neighbors[valid]).min()),
                float((selected[valid] / neighbors[valid]).max()),
            ],
            "eligible_neighbor_count_range": [
                int(neighbors[mask].min()),
                int(neighbors[mask].max()),
            ],
        }

    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "masscubex-mpl"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#aeb8c2",
            "text.color": "#243548",
            "axes.labelcolor": "#243548",
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

    def sticks(ax, mask, color, limit):
        x, y = scan.mz[mask], scan.intensity[mask]
        ax.vlines(x, 0, y, color=color, linewidth=1.3, alpha=0.93, zorder=2)
        above = y > limit
        ax.scatter(
            x[above],
            np.full(above.sum(), limit * 0.975),
            marker="^",
            s=22,
            color=color,
            linewidths=0,
            zorder=5,
        )
        return int(above.sum())

    def panel(ax, row, mask, limit, lo, hi):
        mean, cutoff = (old_mean, old_cutoff) if row == 1 else (new_mean, new_cutoff)
        ax.set(xlim=(lo, hi), ylim=(0, limit), ylabel="Intensity (a.u.)")
        ax.grid(axis="y", color="#e8edf1", linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda value, _: "0" if not value else f"{value / 1000:,.0f}k")
        )
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        if row == 0:
            clipped = sticks(ax, mask, colors["raw"], limit)
        elif row == 1:
            clipped = sticks(ax, mask & old_keep, colors["kept"], limit)
        else:
            clipped = sticks(ax, mask & new_keep & old_keep, colors["kept"], limit)
            clipped += sticks(ax, mask & extra, colors["extra"], limit)
        ax.plot(scan.mz[mask], mean[mask], color=colors["mean"], linewidth=1.5, zorder=3)
        ax.plot(
            scan.mz[mask],
            cutoff[mask],
            color=colors["cutoff"],
            linewidth=1.5,
            linestyle=(0, (4, 2.5)),
            zorder=4,
        )
        hidden = int((mask & (cutoff > limit)).sum())
        label = f"{clipped} tall peaks clipped"
        if hidden:
            label += f"  |  cutoff off-scale at {hidden} points"
        ax.text(
            0.985,
            0.92,
            label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 2},
        )

    handles = [
        Line2D([], [], color=colors["raw"], label="Original centroids"),
        Line2D([], [], color=colors["kept"], label="Kept by comparison rule"),
        Line2D([], [], color=colors["extra"], label="Additional under new rule"),
        Line2D([], [], color=colors["mean"], label="Subset/background mean μ"),
        Line2D(
            [],
            [],
            color=colors["cutoff"],
            linestyle="--",
            label="Active cutoff μ + 3σ",
        ),
    ]
    n = int(window.sum())
    titles = [
        f"Before filtering · {n} centroids (new threshold overlaid)",
        f"{previous_label} · {previous.noise_window_da:g} Da · {int((window & old_keep).sum())}/{n} kept",
        f"Lowest {100 * current.noise_fraction:g}%, min {current.noise_min_neighbors} · {current.noise_window_da:g} Da {current.noise_layout} · {int((window & new_keep).sum())}/{n} kept",
    ]
    fig, axes = plt.subplots(3, 2, figsize=(14, 10.5))
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.17, top=0.80, hspace=0.40, wspace=0.23)
    fig.suptitle(
        (
            f"MS1 background: lowest 5% vs lowest {100 * current.noise_fraction:g}%"
            if args.comparison == "five-percent"
            else f"MS1 background: previous rule vs lowest {100 * current.noise_fraction:g}%"
        ),
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
    for x, heading, value in (
        (0.075, "Comparison cutoff scale", np.median(old_cutoff[window])),
        (0.575, "New cutoff scale", np.median(new_cutoff[window])),
    ):
        fig.text(x, 0.885, heading, fontsize=16, weight="bold")
        fig.text(
            x,
            0.856,
            f"Median = {value / 1000:.1f}k; {100 / args.axis_multiple:.0f}% of axis height",
            fontsize=11,
            color="#526170",
        )
    for row in range(3):
        for col in range(2):
            panel(axes[row, col], row, window, limits[col], args.mz_min, args.mz_max)
            axes[row, col].set_title(titles[row], loc="left", weight="bold", fontsize=12, pad=10)
            if row == 2:
                axes[row, col].set_xlabel("m/z")
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
        "Identical raw intensities in both columns; shared axes within each column. ▲ marks peaks taller than the axis.",
        fontsize=10,
        color="#526170",
    )
    c = counts["shown_window"]
    selected_range = c["selected_count_range_supported"]
    count_label = (
        str(selected_range[0])
        if selected_range[0] == selected_range[1]
        else f"{selected_range[0]}–{selected_range[1]}"
    )
    fig.text(
        0.075,
        0.032,
        f"Full-scan estimates. Current: {current.noise_window_da:g} Da, {count_label} weak neighbors here. Comparison: {previous.noise_window_da:g} Da. Removed ≠ confirmed noise.",
        fontsize=10,
        color="#526170",
    )
    fig.savefig(args.out_dir / "ms1_noise_comparison.png", dpi=160)
    fig.savefig(args.out_dir / "ms1_noise_comparison.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(14, 9))
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.16, top=0.86, hspace=0.42)
    fig.suptitle(
        "Full-scan context · previous and new background rules",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=20,
        weight="bold",
    )
    fig.text(
        0.075,
        0.932,
        f"{args.input.name}  •  scan {scan.native_id}  •  RT {scan.rt:.3f} s",
        fontsize=12,
    )
    for row, (ax, mask, name) in enumerate(
        zip(axes, [all_points, old_keep, new_keep], ["Before", previous_label, "New low-tail"])
    ):
        panel(ax, row, all_points, limits[0], float(scan.mz[0]), float(scan.mz[-1]))
        ax.axvspan(args.mz_min, args.mz_max, color="#dceaf1", alpha=0.45, zorder=0)
        ax.set_title(
            f"{name} · {int(mask.sum())}/{len(scan.mz)} centroids", loc="left", weight="bold", pad=9
        )
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
        "Shaded band = detailed view. Unsupported windows use the zero floor; no noise estimate is inferred from sparse data.",
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
                "previous_mean",
                "previous_sd",
                "previous_cutoff",
                "low_tail_mean",
                "low_tail_sd",
                "low_tail_cutoff",
                "eligible_neighbors",
                "selected_count",
                "kept_previous",
                "kept_low_tail",
            ]
        )
        writer.writerows(
            zip(
                scan.mz,
                scan.intensity,
                old_mean,
                old_sd,
                old_cutoff,
                new_mean,
                new_sd,
                new_cutoff,
                neighbors,
                selected,
                old_keep.astype(int),
                new_keep.astype(int),
            )
        )
    metadata = {
        "version": __version__,
        "source_file": args.input.name,
        "source_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "selection": "nearest acquisition midpoint unless explicit RT; fixed requested m/z range",
        "requested_rt_seconds": requested,
        "acquisition_rt_range_seconds": time_range,
        "native_id": scan.native_id,
        "spectrum_index": scan.spectrum_index,
        "ms1_ordinal_zero_based": ordinal,
        "rt_seconds": scan.rt,
        "centroided": scan.centroided,
        "polarity": scan.polarity,
        "mz_window": [args.mz_min, args.mz_max],
        "config_previous": previous.to_dict(),
        "config_low_tail": current.to_dict(),
        "implementation_sha256": implementation_id()[0],
        "median_previous_cutoff": float(np.median(old_cutoff[window])),
        "median_low_tail_cutoff": float(np.median(new_cutoff[window])),
        "median_low_tail_background": float(np.median(new_mean[window])),
        "axis_multiple": args.axis_multiple,
        "ymax_by_column": limits,
        "counts": counts,
        "interpretation": "Removed/retained centroids have no ground-truth labels. Low-tail sigma is not the standard deviation of the full noise distribution.",
    }
    (args.out_dir / "ms1_noise_example.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"output_directory": str(args.out_dir.resolve()), "counts": counts}, indent=2))


if __name__ == "__main__":
    main()
