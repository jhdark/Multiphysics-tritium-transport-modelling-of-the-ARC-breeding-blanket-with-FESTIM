import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import morethemes as mt
from data_utils import ensure_ornl_data


def signed_percent_formatter(decimals=0):
    def fmt(x, _pos):
        if x == 0:
            return "0%"
        return f"{x:+.{decimals}f}%"

    return FuncFormatter(fmt)


def load_data(values, path_template, column):
    return [
        np.genfromtxt(path_template.format(value=v), delimiter=",", names=True)[column]
        for v in values
    ]


sc_values = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
delta_values = [10, 1, 0.1, 0.01]

# Overlay the OpenFOAM Sc_t sensitivity study on the Sc_t (left) panels
INCLUDE_OPENFOAM = True
openfoam_orange = "#ff7f0e"

# Convert from molecules to milligrams
tritium_mass_mg = 3.01604928 * 1.660539e-27 * 1000000  # u -> kg -> mg
# Extrapolate from 5° sector to full 360° blanket (extensive quantities only)
n_sectors = 72

# (column, sc_file, delta_file, extensive)
# Concentration is intensive, so it is not scaled by the number of sectors
datasets = [
    ("Total_H_volume_1", "total_inventory", "total_inventory", True),
    ("H_advective_flux_surface_2", "outlet_surface_flux", "outlet_surface_flux", True),
    (
        "Average_H_surface_2",
        "outlet_average_concentration",
        "outlet_average_concentration",
        False,
    ),
]

data_sc, data_delta = [], []
for col, sc_file, delta_file, extensive in datasets:
    scale_factor = tritium_mass_mg * (n_sectors if extensive else 1)
    sc = load_data(
        sc_values, f"results/symmBCs/sc_testing/sc_{{value}}/{sc_file}.csv", col
    )
    dt = load_data(
        delta_values,
        f"results/symmBCs/delta_testing/delta_{{value}}/{delta_file}.csv",
        col,
    )
    data_sc.append([v * scale_factor for v in sc])
    data_delta.append([v * scale_factor for v in dt])

# OpenFOAM Sc_t sensitivity. Columns: Sct, Inv (mg), Out flux (mg/s), Out avg conc (mg/m3).
# Rows are normalised to their own Sc_t=0.5 ref, so the absolute units don't matter here.
of_sc_values, data_of_sc = None, [None, None, None]
if INCLUDE_OPENFOAM:
    ensure_ornl_data()
    of = np.genfromtxt(
        "results/OpenFOAM_results_ORNL/OpenFOAM_ARC_tritData.csv",
        delimiter=",",
        skip_header=1,
    )
    of_sc_values = of[:, 0].tolist()
    # inventory, outlet flux, outlet avg. conc. (matches row_configs order)
    # The ORNL concentration column is pre-scaled by the number of sectors;
    # undo it, since concentration is intensive
    data_of_sc = [of[:, 1].tolist(), of[:, 2].tolist(), (of[:, 3] / n_sectors).tolist()]

mt.set_theme("minimal")
plt.rcParams.update({"font.size": 14})

fig, axes = plt.subplots(3, 2, figsize=(10, 8), sharex="col", sharey="row")

festim_blue = (26 / 255, 72 / 255, 72 / 255)

# Per-row config: ylabel, unit, fmt, ref_sc_offset, ref_delta_offset
row_configs = [
    (
        "Inventory",
        "mg",
        ".0f",
        (-32, 10),
        (12, 7),
    ),
    (
        "Outlet flux",
        "mg/s",
        ".2f",
        (8, 10),
        (0, -10),
    ),
    (
        "Outlet avg. conc.",
        "mg/m3",
        ".2f",
        (-13, 12),
        (0, -10),
    ),
]


def rel_diff(ydata, ref_y):
    return [(y / ref_y - 1) * 100 for y in ydata]


def diff_label(v, decimals=1):
    # Values that round to zero are still nonzero -> flag as below the resolution
    if round(v, decimals) == 0:
        return f"<0.{'0' * (decimals - 1)}1%"
    return f"{v:+.{decimals}f}%"


SC_REF_IDX = 1
DELTA_REF_IDX = 2


def plot_sc_panel(
    ax,
    xvals,
    ydata,
    title,
    unit,
    fmt,
    ref_offset,
    left_offset=(0, 10),
    right_offset=(0, 10),
    of_xvals=None,
    of_ydata=None,
):
    ref_idx = SC_REF_IDX
    ref_y = ydata[ref_idx]
    ydiff = rel_diff(ydata, ref_y)

    ax.axhline(0, color=festim_blue, linestyle="--", linewidth=1.5, alpha=0.6, zorder=0)
    if of_ydata is not None:
        of_ydiff = rel_diff(of_ydata, of_ydata[ref_idx])
        # Thick, semi-transparent so a coincident FESTIM line still shows through
        ax.plot(
            of_xvals,
            of_ydiff,
            color=openfoam_orange,
            linewidth=5,
            alpha=0.5,
            solid_capstyle="round",
        )
        # Label the OpenFOAM curve, just above the line
        xc = 1.05
        ax.annotate(
            "OpenFOAM",
            xy=(xc, np.interp(xc, of_xvals, of_ydiff)),
            xytext=(0, 8),
            textcoords="offset points",
            fontsize=12,
            ha="center",
            va="bottom",
            color=openfoam_orange,
        )
    ax.plot(xvals, ydiff, color="black", linewidth=1.8)
    ax.set_xticks([0.3, 0.8, 1.3])
    ax.set_title(title, loc="left", fontsize=20)

    ref_x = xvals[ref_idx]
    # Annotate the first (initial) and last (final) data points
    ax.scatter([ref_x], [0.0], color=festim_blue, zorder=5)
    ax.scatter([xvals[0], xvals[-1]], [ydiff[0], ydiff[-1]], color="black")
    for idx, offset in ((0, left_offset), (-1, right_offset)):
        ax.annotate(
            diff_label(ydiff[idx]),
            xy=(xvals[idx], ydiff[idx]),
            xytext=offset,
            textcoords="offset points",
            fontsize=12,
            ha="center",
            va="top" if offset[1] < 0 else "baseline",
            annotation_clip=False,
        )
    ax.annotate(
        f"Ref \n{ref_y:{fmt}} {unit}",
        xy=(ref_x, 0.0),
        xytext=ref_offset,
        textcoords="offset points",
        fontsize=12,
        ha="left",
        color=festim_blue,
        alpha=0.55,
    )


def plot_delta_panel(
    ax, xvals, ydata, ref_offset, show_minor_xticks=True, first_offset=(0, 10)
):
    ref_idx = DELTA_REF_IDX
    ref_y = ydata[ref_idx]
    ydiff = rel_diff(ydata, ref_y)

    ax.axhline(0, color=festim_blue, linestyle="--", linewidth=1.5, alpha=0.6, zorder=0)
    ax.plot(xvals, ydiff, color="black")
    ax.set_xscale("log")
    ax.set_xticks([1e-2, 1e-1, 1e0, 1e1])
    if not show_minor_xticks:
        ax.tick_params(axis="x", which="minor", bottom=False)

    ref_x = xvals[ref_idx]
    ax.scatter([ref_x], [0.0], color=festim_blue, zorder=5)
    ax.scatter([xvals[0], xvals[-1]], [ydiff[0], ydiff[-1]], color="black")
    ax.annotate(
        diff_label(ydiff[0]),
        xy=(xvals[0], ydiff[0]),
        xytext=first_offset,
        textcoords="offset points",
        fontsize=12,
        ha="center",
        va="top" if first_offset[1] < 0 else "baseline",
    )
    ax.annotate(
        diff_label(ydiff[-1]),
        xy=(xvals[-1], ydiff[-1]),
        xytext=(20, 10),
        textcoords="offset points",
        fontsize=12,
        ha="center",
    )
    ax.annotate(
        "Ref",
        xy=(ref_x, 0.0),
        xytext=ref_offset,
        textcoords="offset points",
        fontsize=12,
        ha="center",
        va="top" if ref_offset[1] < 0 else "baseline",
        color=festim_blue,
        alpha=0.55,
    )


for i, (
    ylabel,
    unit,
    fmt,
    ref_sc_off,
    ref_dt_off,
) in enumerate(row_configs):
    ax_sc, ax_dt = axes[i]
    plot_sc_panel(
        ax_sc,
        sc_values,
        data_sc[i],
        ylabel,
        unit,
        fmt,
        ref_sc_off,
        left_offset=[(0, -12), (0, 12), (0, -12)][i],
        right_offset=[(0, 12), (0, -14), (0, -14)][i],
        of_xvals=of_sc_values,
        of_ydata=data_of_sc[i],
    )
    plot_delta_panel(
        ax_dt,
        delta_values,
        data_delta[i],
        ref_dt_off,
        show_minor_xticks=(i == 2),
        first_offset=(0, -12) if i == 2 else (0, 10),
    )

    # Symmetric y-limits about 0, sized to the larger deviation of the two panels
    row_diffs = rel_diff(data_sc[i], data_sc[i][SC_REF_IDX]) + rel_diff(
        data_delta[i], data_delta[i][DELTA_REF_IDX]
    )
    if data_of_sc[i] is not None:
        row_diffs += rel_diff(data_of_sc[i], data_of_sc[i][SC_REF_IDX])
    lim = 1.2 if i == 1 else max(abs(min(row_diffs)), abs(max(row_diffs))) * 1.2
    ax_sc.set_ylim(-lim, lim)
    ax_sc.yaxis.set_major_formatter(signed_percent_formatter(0))
    if i == 1:
        ax_sc.set_yticks([-1, 0, 1])
    for label in (
        ax_sc.get_yticklabels() + ax_sc.get_xticklabels() + ax_dt.get_xticklabels()
    ):
        label.set_alpha(0.6)

# Add x-labels only on bottom row
axes[2][0].set_xlabel(r"Turbulent Schmidt No., Sc$_\mathrm{t}$", fontsize=18)
axes[2][1].set_xlabel(r"Tuning parameter, $\delta$", fontsize=18)

plt.tight_layout(h_pad=2.5)

plt.savefig("figures/plot_parametric.pdf", dpi=300, bbox_inches="tight")
plt.savefig("figures/plot_parametric.png")
plt.show()
