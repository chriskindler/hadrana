# /hadrana/src/hadrana/misc/lhpc/lhpc_visualise.py
import matplotlib.pyplot as plt
from cycler import cycler

plt.rcParams.update({
    "font.size": 14,
    "font.family": "computer modern",
    "text.usetex": True,
    "figure.dpi": 200,
    "axes.prop_cycle": cycler(color=[
        "#0091ff", # blue
        "#ff8d28", # orange
        "#34C759", # green
        "#ff2d55", # pink
    ])
})

t_color    = {8: "#ff8d28", 10: "#0091ff", 12: "#ff2d55", 14: "#34C759", 19: "#ff8d28"}
# t_color    = {10: "tab:blue", 12: "tab:orange", 14: "tab:green", 19: "tab:red"}
t_marker   = {8: "d", 10: "o", 12: "^", 14: "v", 19: "o"}
markersize = 6
capsize    = 6
elinewidth = 1.3