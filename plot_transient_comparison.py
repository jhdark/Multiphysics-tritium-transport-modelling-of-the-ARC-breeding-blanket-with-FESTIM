import matplotlib.pyplot as plt
import numpy as np
import morethemes as mt
from data_utils import ensure_ornl_data

# ---------------------------------------------------------------------------
# FESTIM transient data
# ---------------------------------------------------------------------------
surface_flux_data = np.genfromtxt(
    "results/symmBCs/transient/outlet_surface_flux.csv", delimiter=",", names=True
)
festim_flux_t = np.insert(surface_flux_data["ts"] / 3600, 0, 0)  # Convert to hours

total_inventory_data = np.genfromtxt(
    "results/symmBCs/transient/total_inventory.csv", delimiter=",", names=True
)
festim_inv_t = np.insert(total_inventory_data["ts"] / 3600, 0, 0)  # Convert to hours

# Convert from molecules to milligrams and extrapolate from 5° sector to full 360°
scale_factor = 3.01604928 * 1.660539e-27 * 1000000 * 72  # u -> kg -> mg
festim_flux = np.insert(
    surface_flux_data["H_advective_flux_surface_2"] * scale_factor, 0, 0
)
festim_inv = np.insert(total_inventory_data["Total_H_volume_1"] * scale_factor, 0, 0)

# ---------------------------------------------------------------------------
# OpenFOAM transient data (ORNL)
# ---------------------------------------------------------------------------
ensure_ornl_data()
of_path ="results/OpenFOAM_results_ORNL/postProcessing"
# Outlet flux: time | sum(tritFlux)
of_flux_data = np.loadtxt(
    f"{of_path}/tritFlux_outlet/0/surfaceFieldValue.dat", comments="#"
)
# Volume integral: time | volIntegrate(tritConc) | volIntegrate(tritProd)
of_inv_data = np.loadtxt(f"{of_path}/tritVolInt/0/volFieldValue.dat", comments="#")

# 72: 5° sector -> full 360° blanket, 1e6: kg -> mg
of_scale = 72 * 1e6
of_flux_t = of_flux_data[:, 0] / 3600  # Convert to hours
of_flux = of_flux_data[:, 1] * of_scale
of_inv_t = of_inv_data[:, 0] / 3600  # Convert to hours
of_inv = of_inv_data[:, 1] * of_scale

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
festim_blue = (26 / 255, 72 / 255, 72 / 255)
openfoam_orange = "#ff7f0e"

mt.set_theme("minimal")
plt.rcParams.update({"font.size": 14})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), sharex=True)

x_label = 0.5  # hours, where the curve-name labels are centered
t_max = 1.0  # hours (matches ax xlim); end markers use last point within this window


def end_index(t):
    return np.searchsorted(t, t_max, side="right") - 1


def label_curve(ax, t, y, name, color, above):
    """Annotate a curve with its name near x=x_label, above or below the line."""
    dy = 6 if above else -6
    va = "bottom" if above else "top"
    ax.annotate(
        name,
        xy=(x_label, np.interp(x_label, t, y)),
        xytext=(0, dy),
        textcoords="offset points",
        fontsize=14,
        ha="center",
        va=va,
        color=color,
    )


def label_end(ax, t, y, text, color, above):
    """Dot + final-value annotation at the end of the curve, above or below it."""
    i = end_index(t)
    x_end = 0.99
    ax.scatter([x_end], [y[i]], color=color, zorder=5)
    dy = 8 if above else -8
    va = "bottom" if above else "top"
    ax.annotate(
        text,
        xy=(x_end, y[i]),
        xytext=(-6, dy),
        textcoords="offset points",
        fontsize=12,
        ha="right",
        va=va,
        color=color,
    )


# Top subplot: Outlet flux (OpenFOAM above, FESTIM below)
ax1.plot(of_flux_t, of_flux, color=openfoam_orange)
ax1.plot(festim_flux_t, festim_flux, color=festim_blue)
ax1.set_ylim(0, 1.2)
ax1.set_yticks([0, 0.4, 0.8, 1.2])
ax1.set_title(
    r"Outlet flux (mg$\,\mathrm{s}^{-1}$)", fontsize=18, loc="left", color="black"
)
label_curve(ax1, of_flux_t, of_flux, "OpenFOAM", openfoam_orange, above=True)
label_curve(ax1, festim_flux_t, festim_flux, "FESTIM", festim_blue, above=False)
label_end(
    ax1,
    of_flux_t,
    of_flux,
    rf"{of_flux[end_index(of_flux_t)]:.2f} mg$\,\mathrm{{s}}^{{-1}}$",
    openfoam_orange,
    above=True,
)
label_end(
    ax1,
    festim_flux_t,
    festim_flux,
    rf"{festim_flux[end_index(festim_flux_t)]:.2f} mg$\,\mathrm{{s}}^{{-1}}$",
    festim_blue,
    above=False,
)

# Bottom subplot: Total inventory (FESTIM above, OpenFOAM below)
ax2.plot(of_inv_t, of_inv, color=openfoam_orange)
ax2.plot(festim_inv_t, festim_inv, color=festim_blue)
ax2.set_xlabel("Time (h)", fontsize=18, color="black")
ax2.set_ylim(0, 300)
ax2.set_xlim(0, 1)
ax2.set_xticks([0, 0.5, 1.0])
ax2.set_yticks([0, 100, 200, 300])
ax2.set_title("Inventory (mg)", fontsize=18, loc="left", color="black")
label_curve(ax2, festim_inv_t, festim_inv, "FESTIM", festim_blue, above=True)
label_curve(ax2, of_inv_t, of_inv, "OpenFOAM", openfoam_orange, above=False)
label_end(
    ax2,
    festim_inv_t,
    festim_inv,
    f"{festim_inv[end_index(festim_inv_t)]:.0f} mg",
    festim_blue,
    above=True,
)
label_end(
    ax2,
    of_inv_t,
    of_inv,
    f"{of_inv[end_index(of_inv_t)]:.0f} mg",
    openfoam_orange,
    above=False,
)

# Show a horizontal gridline at the bottom (y=0) but not at the top tick
for ax in (ax1, ax2):
    gridlines = ax.get_ygridlines()
    gridlines[0].set_clip_on(False)  # bottom (0) line: render fully at the edge
    gridlines[-1].set_visible(False)  # hide the top gridline

plt.tight_layout()
plt.savefig("results/plot_transient_comparison.png")
plt.savefig("results/plot_transient_comparison.pdf")

plt.show()
