from cycler import cycler
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

from cora.ensembles.helpers import *
from cora.momenta import *

ensemble = "H107"
Nbst = 100

def load_ratios(ensemble, source_set, kind):
    path = f"/glurch/scratch/kic04594/data/misc/{ensemble}-{kind}-ISOVECTOR-AXIAL-RATIOS-NBST{Nbst}-{source_set}.pkl"
    with open(path, "rb") as f:
        return pickle.load(f)

ratios_spat_s2 = load_ratios(ensemble, "source_set2", "SPATIAL")
ratios_spat_s4 = load_ratios(ensemble, "source_set4", "SPATIAL")

qmax = 8
qcol: dict[str] = momentum_index_dict(qmax)              # e.g. "(0, 0, 1)" -> idx
qidx_to_qvec = {idx: qstr for qstr, idx in qcol.items()} # idx -> "(0, 0, 1)"
qidx_list = list(qcol.values())

plt.rcParams.update({
    "font.size": 14,
    "font.family": "computer modern",
    "text.usetex": True,
    "figure.dpi": 200,
    "axes.prop_cycle": cycler(color=[
        "#0091ff", "#ff8d28", "#34C759", "#ff2d55",
    ])
})

ipol_X = 0
jcur_x = 0
ts = [10, 12, 14]
src_ids_to_plot = [0, 1]

def mean_err_over_bootstrap(arr):
    mu = arr.mean(axis=0)
    err = arr.std(axis=0, ddof=1)
    return mu, err

def get_tau_curve(ratios_dict, t, source_set, src_id, ipol, jcur, qidx):
    A = ratios_dict[t][source_set]  # (Nbst, Npol, Ncur, Nsrc, Nmom, Ntau)
    return A[:, ipol, jcur, src_id, qidx, :]  # (Nbst, Ntau)

def tau_axis(Ntau, centered: bool):
    if centered:
        return np.arange(Ntau) - (Ntau - 1) / 2
    return np.arange(Ntau)

t_color  = {10: "#0091ff", 12: "#ff2d55", 14: "#34C759"}
t_marker = {10: "o", 12: "^", 14: "v"}
markersize = 6
capsize = 6
elinewidth = 1.3

def sanitize_qvec(qvec_str: str) -> str:
    # "(0, 0, 1)" -> "0_0_1"
    return qvec_str.strip().replace("(", "").replace(")", "").replace(" ", "").replace(",", "_")

outdir = f"/glurch/home/kic04594/hadrana/src/hadrana/misc/lhpc/plots/{ensemble}/POLX_CURX"
os.makedirs(outdir, exist_ok=True)

for qidx in qidx_list:
    qvec = qidx_to_qvec[qidx]
    qtag = sanitize_qvec(qvec)

    # NEW FIGURE PER MOMENTUM
    fig, axes = plt.subplots(
        nrows=len(src_ids_to_plot) * len(ts),  # 6 rows
        ncols=2,
        figsize=(11, 18),
        sharex=False,
        sharey=True,
        constrained_layout=True,
    )

    for it, t in enumerate(ts):
        colour = t_color[t]

        for row, src_id in enumerate(src_ids_to_plot):
            global_row = it * len(src_ids_to_plot) + row

            for col, (label, ratios, source_set) in enumerate([
                ("LHP", ratios_spat_s2, "source_set2"),
                ("REG", ratios_spat_s4, "source_set4"),
            ]):
                ax = axes[global_row, col]

                # guards (optional but safer)
                if t not in ratios or source_set not in ratios[t]:
                    ax.axis("off")
                    continue
                if src_id >= ratios[t][source_set].shape[3]:
                    ax.axis("off")
                    continue
                if qidx >= ratios[t][source_set].shape[4]:
                    ax.axis("off")
                    continue

                y_bst = get_tau_curve(ratios, t, source_set, src_id, ipol_X, jcur_x, qidx)
                est, err = mean_err_over_bootstrap(y_bst)
                tau = tau_axis(len(est), centered=(qidx == 0))

                ax.errorbar(
                    tau[1:-1],
                    est[1:-1],
                    yerr=err[1:-1],
                    marker=t_marker[t],
                    markersize=markersize,
                    linestyle=":",
                    color=colour,
                    elinewidth=elinewidth,
                    capsize=capsize,
                    label=rf"$t={t}a$",
                )

                ax.grid(True, linestyle=":", linewidth=0.5)

                ax.set_title(
                    rf"$\mathsf{{{ensemble}}},$ {label}, source id = {src_id}, "
                    rf"$\mathbf{{q}} \sim {qvec}^{{\mathsf{{T}}}}$"
                )

                if global_row == len(ts) * len(src_ids_to_plot) - 1:
                    ax.set_xlabel(r"$\tau/a - t/2a$" if qidx == 0 else r"$\tau/a$")

                if col == 0:
                    ax.set_ylabel(r"$\mathrm{Im}\,R_A^{x,x}(\mathbf{q}, t, \tau)$")

                ax.legend()

    filepath = f"{outdir}/{ensemble}-ISOVECTOR-AXIAL-RATIOS-LHP-REG-POLX-CURX-MOM-{qtag}.pdf"
    print(f"[SAVE] qidx={qidx} qvec={qvec} -> {filepath}")
    fig.savefig(filepath, bbox_inches="tight")
    plt.close(fig)