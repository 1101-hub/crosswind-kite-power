"""
The design chart: given a kite area and a wind regime, what do you get?

This is the deliverable the whole model exists to produce. Two figures:

  1. DIMENSIONAL  -- power against kite area and wind speed, with the
     no-flight boundary drawn on it, and real systems marked. Spans from a
     festival patang to the Makani M600, a factor of 450 in area and roughly
     250,000 in power.

  2. NON-DIMENSIONAL -- the same physics collapsed onto two numbers, the
     glide ratio E and the gravity number G. This is the version that
     generalises: any kite with the same (E, G) sits at the same point,
     regardless of size.

The no-flight boundary is the interesting feature. It comes from the
discriminant of the flight-speed quadratic going negative, and it is what
separates "a small kite in light wind" from "a small kite that cannot fly
crosswind at all".

Colour follows the dataviz rules: magnitude gets ONE hue, light to dark. No
rainbow -- a rainbow ramp invents boundaries in the data that are not there.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

import kite as K


# Single-hue sequential ramp, matching the interactive page.
RAMP = LinearSegmentedColormap.from_list(
    "tension", ["#FDF0F5", "#FBDCE6", "#F7A8C2", "#F0407A", "#B41E52", "#741132"])

INK, INK2, INK3 = "#12181C", "#46545C", "#74848C"
LINE, FLOW = "#DDE2E4", "#00A0A8"

REAL_SYSTEMS = [
    ("Patang",  0.12, 0.012, 0.70, 0.25, 0.0005, 1.6,  50.0),
    ("Trainer", 1.5,  0.9,   0.90, 0.20, 0.0015, 4.0,  50.0),
    ("Foil 3m²", 3.0, 1.8,   0.90, 0.20, 0.002,  4.0, 100.0),
    ("Makani M600", 54.0, 1850.0, 1.4, 0.10, 0.03, 12.0, 440.0),
]


def make_cfg(area, mass, cl, cd, d, ar, L):
    return {'area': area, 'mass': mass, 'cl': cl, 'cd': cd, 'cl_depowered': 0.15,
            'tether_diameter': d, 'aspect_ratio': ar, 'tether_length': L,
            'zeta_turn': 1.0, 'tension_limit': 1e9}


def best_power(cfg, v, n=96, want_pattern=False):
    """Coarse optimisation over BOTH path families at one operating point.

    A full Nelder-Mead run per grid cell would take hours for a whole chart.
    The optimum is smooth and stable, so a modest grid is enough -- and a grid
    cannot get stuck in a local optimum, which matters when sweeping across
    the stall boundary.

    Both a figure-eight and a circle are tried, because which one is even
    FEASIBLE depends on size. The five-wingspan turning limit scales as
    sqrt(area), while the tightest turn a Lissajous eight can hold is about
    0.18 of the tether length against a circle's 0.22. Past a certain scale
    the eight is geometrically impossible and only the circle survives.
    """
    span = np.sqrt(cfg['aspect_ratio'] * cfg['area'])
    r_min = 5.0 * span
    best, winner = 0.0, "none"

    reels = (0.10, 0.16, 0.22, 0.28, 0.34)

    for beta0 in (18.0, 25.0, 32.0, 40.0, 50.0):
        # --- figure-eight ------------------------------------------------
        for sweep in (20.0, 32.0, 44.0, 56.0):
            for rise in (6.0, 12.0, 18.0):
                if rise >= beta0 - 13.0:
                    continue
                path = {'beta0': beta0, 'sweep': sweep, 'rise': rise}
                for f in reels:
                    p, _, _, stall, rturn = K.reel_out_power(
                        v, cfg, path, f, n=n, state_fn=K.lissajous_state)
                    if stall > 0.02 or rturn < r_min * 0.98:
                        continue
                    if p > best:
                        best, winner = p, "eight"

        # --- circle -------------------------------------------------------
        for R in (10.0, 15.0, 20.0, 25.0, 30.0):
            if R >= beta0 - 8.0:
                continue
            path = {'beta0': beta0, 'radius': R}
            for f in reels:
                p, _, _, stall, rturn = K.reel_out_power(
                    v, cfg, path, f, n=n, state_fn=K.circle_state)
                if stall > 0.02 or rturn < r_min * 0.98:
                    continue
                if p > best:
                    best, winner = p, "circle"

    return (best, winner) if want_pattern else best


def dimensional_chart(path_out="../figures/design_chart.png"):
    areas = np.logspace(np.log10(0.05), np.log10(80), 26)
    winds = np.linspace(1.5, 14.0, 24)

    P = np.zeros((len(winds), len(areas)))
    W = np.zeros((len(winds), len(areas)))       # 1 = eight, 2 = circle

    for i, v in enumerate(winds):
        for j, A in enumerate(areas):
            # Scale the kite plausibly with area: soft-kite areal density,
            # tether diameter set by the load it carries, and tether LENGTH
            # set by wingspan.
            #
            # That last one matters more than it looks. The turning limit is
            # 5 wingspans and the achievable turn radius is a fixed fraction
            # of tether length, so a tether that does not grow with span makes
            # large kites geometrically unflyable. An earlier version scaled
            # L as A^0.16 and produced a chart where nothing above 1 m2 could
            # fly -- which was the model correctly reporting an impossible
            # configuration, not a bug in the physics.
            span = np.sqrt(3.5 * A)
            mass = 0.10 * A ** 1.15
            d = 0.0006 * np.sqrt(A / 0.12)
            L = max(30.0, 25.0 * span)
            cfg = make_cfg(A, mass, 0.85, 0.22, d, 3.5, L)
            p, pat = best_power(cfg, v, want_pattern=True)
            P[i, j] = p
            W[i, j] = {"eight": 1.0, "circle": 2.0}.get(pat, 0.0)

    fig, ax = plt.subplots(figsize=(9.2, 6.0), dpi=200)
    fig.patch.set_facecolor("white")

    masked = np.ma.masked_where(P <= 1e-4, P)
    lo = max(masked.min(), 1e-3)
    cs = ax.contourf(areas, winds, masked, levels=np.logspace(
        np.log10(lo), np.log10(masked.max()), 22), norm=LogNorm(), cmap=RAMP)

    cl = ax.contour(areas, winds, masked, levels=[1, 10, 100, 1e3, 1e4, 1e5],
                    colors=INK, linewidths=.6, alpha=.45)
    ax.clabel(cl, fmt=lambda x: ("%g W" % x) if x < 1000 else ("%g kW" % (x / 1000)),
              fontsize=7.5, inline=True)

    # the no-flight region
    ax.contourf(areas, winds, P, levels=[-1, 1e-4], colors=["#EFF2F3"])
    ax.contour(areas, winds, P, levels=[1e-4], colors=[FLOW], linewidths=1.6)
    # Label the no-flight region wherever it actually is. It moves: small
    # kites are light and fly in very light wind, while large kites carry
    # enough mass that gravity beats them at the bottom of the wind range.
    if (P <= 1e-4).any():
        ii, jj = np.where(P <= 1e-4)
        k = int(np.argmax(areas[jj]))
        ax.annotate("no sustained crosswind flight\n(quadratic has no real root)",
                    xy=(areas[jj[k]], winds[ii[k]]),
                    xytext=(-10, 34), textcoords="offset points",
                    fontsize=8, color=FLOW, style="italic", ha="right",
                    arrowprops=dict(arrowstyle="-", color=FLOW, lw=.9))

    for name, A, m, cl_, cd_, d, ar, L in REAL_SYSTEMS:
        if A < areas[0] or A > areas[-1]:
            continue
        ax.plot([A], [6.0], "o", ms=6, mfc="white", mec=INK, mew=1.4, zorder=6)
        ax.annotate(name, (A, 6.0), textcoords="offset points", xytext=(0, 11),
                    ha="center", fontsize=8.5, color=INK, weight="600", zorder=6)

    ax.set_xscale("log")
    ax.set_xlabel("Kite area  (m²)", fontsize=10, color=INK2, labelpad=8)
    ax.set_ylabel("Wind speed  (m/s)", fontsize=10, color=INK2, labelpad=8)
    ax.set_title("Cycle-averaged power from a crosswind kite",
                 fontsize=13.5, color=INK, weight="600", loc="left", pad=30)
    ax.text(0, 1.012, "best of figure-eight and circular loop at each point · quasi-steady model · "
            "markers placed at 6 m/s for reference",
            transform=ax.transAxes, fontsize=8.5, color=INK3)

    ticks = [t for t in (0.1, 1, 10, 100, 1e3, 1e4, 1e5)
             if masked.min() <= t <= masked.max()]
    cb = fig.colorbar(cs, ax=ax, pad=.02, ticks=ticks)
    cb.ax.set_yticklabels([("%g W" % t) if t < 1000 else ("%g kW" % (t / 1000))
                           for t in ticks])
    cb.set_label("cycle-averaged power", fontsize=9, color=INK2)
    cb.outline.set_edgecolor(LINE)

    for sp in ax.spines.values():
        sp.set_color(LINE)
    ax.tick_params(colors=INK3, labelsize=8.5)
    fig.tight_layout()
    fig.savefig(path_out, facecolor="white", bbox_inches="tight")
    print("wrote", path_out)
    return areas, winds, P


def nondimensional_chart(path_out="../figures/design_chart_nondim.png"):
    """Power coefficient against the two numbers that actually govern it.

        E = C_L / C_D_total          how good a wing it is
        G = mg / (0.5 rho v^2 A C_L) how heavy it is, relative to its lift

    Any two kites sharing (E, G) behave identically here, whatever their size.
    That is the whole point of non-dimensionalising, and it is what lets a
    0.12 m² patang say anything about a 54 m² Makani wing.
    """
    Es = np.linspace(1.5, 12.0, 26)
    Gs = np.logspace(np.log10(0.005), np.log10(1.2), 24)

    Z = np.zeros((len(Gs), len(Es)))
    A, v, cl = 3.0, 6.0, 0.9
    q = 0.5 * K.RHO * v ** 2 * A * cl

    for i, Gn in enumerate(Gs):
        for j, E in enumerate(Es):
            mass = Gn * q / K.G
            cd_tot = cl / E
            cfg = make_cfg(A, mass, cl, cd_tot, 1e-9, 4.0, 100.0)
            p = best_power(cfg, v)
            Z[i, j] = p / (0.5 * K.RHO * A * v ** 3)

    fig, ax = plt.subplots(figsize=(8.4, 5.6), dpi=200)
    fig.patch.set_facecolor("white")

    masked = np.ma.masked_where(Z <= 1e-5, Z)
    cs = ax.contourf(Es, Gs, masked, levels=24, cmap=RAMP)
    cl2 = ax.contour(Es, Gs, masked, levels=8, colors=INK, linewidths=.6, alpha=.45)
    ax.clabel(cl2, fontsize=7.5, fmt="%.2f", inline=True)

    ax.contourf(Es, Gs, Z, levels=[-1, 1e-5], colors=["#EFF2F3"])
    ax.contour(Es, Gs, Z, levels=[1e-5], colors=[FLOW], linewidths=1.6)

    # The 1.5 m2 trainer and the 3 m2 foil share (E, G) exactly, so they land
    # on the same point despite one being twice the area. That is the collapse
    # doing its job -- label them together rather than stacking two markers.
    for name, Ev, Gv in (("Patang", 2.3, 0.064),
                         ("Trainer & 3 m² foil  (same E, G)", 4.2, 0.30),
                         ("Makani M600", 11.0, 0.35)):
        if Es[0] <= Ev <= Es[-1] and Gs[0] <= Gv <= Gs[-1]:
            ax.plot([Ev], [Gv], "o", ms=6, mfc="white", mec=INK, mew=1.4, zorder=6)
            ax.annotate(name, (Ev, Gv), textcoords="offset points", xytext=(7, 5),
                        fontsize=8.5, color=INK, weight="600", zorder=6)

    ax.set_yscale("log")
    ax.set_xlabel("Glide ratio  E = C_L / C_D", fontsize=10, color=INK2, labelpad=8)
    ax.set_ylabel("Gravity number  G = mg / (½ρv²A·C_L)", fontsize=10, color=INK2, labelpad=8)
    ax.set_title("Power coefficient  P / (½ρAv³)", fontsize=13.5, color=INK,
                 weight="600", loc="left", pad=30)
    ax.text(0, 1.012, "the same physics, size removed — any kite sharing (E, G) "
            "lands on the same point", transform=ax.transAxes, fontsize=8.5, color=INK3)

    cb = fig.colorbar(cs, ax=ax, pad=.02)
    cb.set_label("P / (½ρAv³)", fontsize=9, color=INK2)
    cb.outline.set_edgecolor(LINE)

    for s in ax.spines.values():
        s.set_color(LINE)
    ax.tick_params(colors=INK3, labelsize=8.5)
    fig.tight_layout()
    fig.savefig(path_out, facecolor="white", bbox_inches="tight")
    print("wrote", path_out)


if __name__ == "__main__":
    import os
    os.makedirs("../figures", exist_ok=True)
    dimensional_chart()
    nondimensional_chart()
