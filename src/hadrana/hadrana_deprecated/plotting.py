import gvar as gv
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FixedLocator, FormatStrFormatter, FuncFormatter

from hadrana.ensembles import *
from hadrana.serialiser import *

def plot_form_factor_data(
        ensemble_ids: list[str],
        data,
    ):

    ids_sym = ["X250", "X251", "D250"]
    ids_trm = ["N203", "D251", "D200"]
    ids_msc = ["N204", "N201", "D201"]

    marker_map = {
        # HEAVY 
        "X250": "s",
        "N203": "D",
        "N204": "8",
        # MID 
        "X251": "v",
        "D251": "^",
        "N201": "<",
        # LIGHT
        "D250": "d",
        "D200": "h",
        "D201": "p"
    }

    try:
        marker_map
    except NameError:
        marker_map = {
            "X250": "o", "X251": "s", "D250": "D",
            "N203": "^", "D251": "v", "D200": "P",
            "N204": "X", "N201": "*", "D201": "h",
        }

    Q2      = np.stack([data[eid]["Q"] for eid in ensemble_ids], axis=0)
    GA_gv   = [data[eid]["GA"] for eid in ensemble_ids]
    GA_cen  = np.stack([gv.mean(arr) for arr in GA_gv], axis=0)
    GA_err  = np.stack([gv.sdev(arr) for arr in GA_gv], axis=0)

    map_sym = {"X250": 348, "X251": 267, "D250": 197}
    map_trm = {"N203": 346, "D251": 286, "D200": 201}
    map_msc = {"N204": 351, "N201": 285, "D201": 199}
    all_maps = {"SYM": map_sym, "TRM": map_trm, "MSC": map_msc}
    color_map = {"SYM": "tab:blue", "TRM": "tab:orange", "MSC": "tab:green"}

    fig, axs = plt.subplots(1, 2, figsize=(9, 4), gridspec_kw={'width_ratios': [3, 1]})

    for i, eid in enumerate(ensemble_ids):
        if eid in ids_sym:
            color = "tab:blue"
        elif eid in ids_trm:
            color = "tab:orange"
        elif eid in ids_msc:
            color = "tab:green"
        else:
            color = "tab:gray"

        axs[0].errorbar(
            Q2[i], GA_cen[i], yerr=GA_err[i],
            markersize=7, capsize=6, color=color, fmt=marker_map[eid],
            label=eid if (Q2.shape[1] and i == 0) else None
        )

    axs[0].grid(True, which='both', alpha=1)
    axs[0].set_xlabel(r"$Q^2$ [GeV$^2$]")
    axs[0].set_ylabel(r"$G_A(Q^2)$")


    all_mps = list(map_sym.values()) + list(map_trm.values()) + list(map_msc.values())
    ymin, ymax = min(all_mps), max(all_mps)

    for group, group_map in all_maps.items():
        x_coord = r"$\mathsf{" + group + r"}$"
        for ens_id, mp in group_map.items():
            axs[1].scatter(x=x_coord, y=mp, color=color_map[group],
                        marker=marker_map.get(ens_id, "o"), s=80)
            axs[1].text(x=x_coord, y=mp + 6, s=rf"$\mathsf{{{ens_id}}}$",
                        ha='center', va='bottom', fontsize=12)

    axs[0].grid(True, alpha=0.3)

    for ax in axs:
        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)

    axs[1].set_ylabel(r"$M_{\pi}$ [MeV]", rotation=270, labelpad=20)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.tick_right()
    axs[1].set_ylim(ymin - 25, ymax + 25)
    axs[1].margins(x=0.2)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.10)  
    plt.show()