import shutil
import matplotlib.pyplot as plt
from cycler import cycler

USE_TEX = shutil.which("latex") is not None

plt.rcParams.update({
    "font.size": 15,
    "figure.dpi": 400,
    "axes.prop_cycle": cycler(color=["tab:blue","tab:orange","tab:green","tab:red"]),
    "text.usetex": USE_TEX,
})

if USE_TEX:
    plt.rcParams.update({
        "font.family": "serif",
    })
else:
    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "cm",
    })