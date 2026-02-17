# DEFAULT
import matplotlib.pyplot as plt
from cycler import cycler

plt.rcParams.update({
    "font.size": 14,
    "font.family": "computer modern",
    "text.usetex": True,
    "figure.dpi": 400,

    "axes.prop_cycle":
        cycler(color=[
            # default
            "tab:blue",
            "tab:orange",
            "tab:green",
            "tab:red",
            # apple dev colors:
            # "#0088ff",
            # "#ff8d28",
            # "#34c759",
            # "#ff383c",

             "#1e6ef4", # blue:   Increased contrast (light)
            "#e7124d", # pink:   Increased contrast (light)
            "#008575", # mint:   Increased contrast (light)
            "#ff8d28", # orange: Default (light)
        ])
    })


