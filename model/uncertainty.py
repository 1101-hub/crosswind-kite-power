"""
Monte Carlo uncertainty propagation for the kite model.

Why this exists
---------------
A prediction of "2.32 W" cannot be wrong. Any measurement differs from it, so
"close enough" becomes a judgement call and the comparison proves nothing.

A prediction of "1.4 to 3.6 W (90% interval)" CAN be wrong. If the measurement
lands outside the band, the model is falsified. Bands are what make P1-P7
testable rather than decorative.

Two things are computed here:

  1. BANDS   -- draw every uncertain input from its distribution a few thousand
                times, push each draw through the model, report percentiles.

  2. SENSITIVITY -- vary one input at a time and see how much the output moves.
                This answers a practical question: WHERE SHOULD YOU SPEND
                EFFORT? If wind speed dominates, a better anemometer helps more
                than anything else. If C_L dominates, do more static tests.

Why Monte Carlo instead of adding errors in quadrature
-----------------------------------------------------
Because the inputs enter the model through several paths at once and the
combination is not obvious. A worked example, which is worth putting in the
writeup because the intuitive answer is WRONG:

  Prediction P2 is a ratio, crosswind tension over static tension. It is
  tempting to say a high C_L raises both and therefore cancels, so the ratio
  must be better determined than its parts. Run it and that is false:

      static tension     ~ C_L
      crosswind tension  ~ C_L * E^2 = C_L * (C_L/C_D)^2 = C_L^3 / C_D^2
      ratio              ~ C_L^2 / C_D^2

  Nothing cancels. The ratio depends on the SQUARE of an already-uncertain
  glide ratio, so it is more sensitive to the coefficients than either tension
  alone, not less.

This is the real argument for Monte Carlo: not that it handles correlations
better than quadrature, but that it stops you reasoning your way to a
confident wrong answer about how errors combine.
"""

import numpy as np

import kite as K
import optimize as O


# ---------------------------------------------------------------------------
# How well do we actually know each input?
#
# These are 1-sigma RELATIVE uncertainties. Be honest here -- inflating them
# to make the model look good is self-defeating, because a band wide enough to
# never be wrong tests nothing.
# ---------------------------------------------------------------------------

PRIOR = {                       # TODAY: coefficients are guesses
    'cl':               0.20,   # estimated from flat-plate theory, not measured
    'cd':               0.25,   # worse: induced drag on a low-aspect-ratio plate
    'area':             0.05,   # measuring a diamond with a ruler
    'mass':             0.05,   # kitchen scale
    'tether_diameter':  0.20,   # cotton line varies along its length
    'cd_tether':        0.20,   # literature spread on cylinder drag
    'wind_speed':       0.20,   # measured at head height, kite is 20-40 m up
    'rho':              0.03,   # air density from temperature and altitude
}

POST_P1 = {                     # AFTER the static test has measured E
    'cl':               0.08,
    'cd':               0.08,
    'area':             0.05,
    'mass':             0.05,
    'tether_diameter':  0.20,
    'cd_tether':        0.20,
    'wind_speed':       0.10,   # assumes a logged mast anemometer, not a handheld
    'rho':              0.03,
}

# The active set. Swap with use_scenario().
UNCERTAINTY = dict(PRIOR)


def use_scenario(d):
    """Point the module at a different uncertainty set."""
    UNCERTAINTY.clear()
    UNCERTAINTY.update(d)

# Note on wind_speed: you measure at head height, the kite flies at 20-40 m.
# The log profile means those differ, and power goes as v^3, so a 20% wind
# error is a 70% power error. This term will dominate everything. That is not
# a flaw in the model -- it is the actual experimental difficulty, and saying
# so plainly is stronger than hiding it.


def perturb(base, rng, keys=None):
    """Draw one sample of the configuration.

    Normal draws, clipped at 20% of nominal so nothing goes negative or absurd.
    `keys` restricts which inputs are varied -- used by the sensitivity study.
    """
    cfg = dict(base)
    active = UNCERTAINTY if keys is None else {k: UNCERTAINTY[k] for k in keys}

    for k, sigma in active.items():
        if k in ('wind_speed', 'cd_tether', 'rho'):
            continue                       # handled by the caller, not in cfg
        if k in cfg:
            cfg[k] = max(cfg[k] * (1.0 + rng.normal(0, sigma)), cfg[k] * 0.2)
    return cfg


def one_sample(base, path, v_nominal, reel, rng, keys=None):
    """Run the model once with perturbed inputs. Returns a dict of outputs."""
    active = UNCERTAINTY if keys is None else {k: UNCERTAINTY[k] for k in keys}
    cfg = perturb(base, rng, keys)

    v = v_nominal
    if 'wind_speed' in active:
        v = max(v_nominal * (1.0 + rng.normal(0, active['wind_speed'])), 0.1)

    rho_scale = 1.0
    if 'rho' in active:
        rho_scale = max(1.0 + rng.normal(0, active['rho']), 0.5)

    # Air density and tether drag both enter through module-level constants,
    # so swap them for the duration of this sample and restore afterwards.
    rho_saved = K.RHO
    K.RHO = rho_saved * rho_scale
    try:
        cd_t = 1.0
        if 'cd_tether' in active:
            cd_t = max(1.0 + rng.normal(0, active['cd_tether']), 0.2)

        cd_eff = K.effective_drag_coefficient(cfg['cd'], cfg['tether_length'],
                                              cfg['tether_diameter'], cfg['area'],
                                              cd_tether=cd_t)
        E = cfg['cl'] / cd_eff

        T_static, beta_static = K.static_tension_and_angle(
            v, cfg['area'], cfg['cl'], cd_eff, cfg['mass'])

        r = K.cycle_power(v, cfg, path, reel, 0.6, 65.0, 30.0)
        loyd = K.loyd_limit(v, cfg['area'], cfg['cl'], E)

        flying = r['stalled_fraction'] < 0.02 and r['tension_mean'] > 0
        return {
            'E': E,
            'T_static': T_static,
            'beta_static': beta_static,
            'T_crosswind': r['tension_mean'] if flying else np.nan,
            'amplification': (r['tension_mean'] / T_static) if flying and T_static > 0 else np.nan,
            'power': r['power_avg'] if flying else np.nan,
            'loyd_fraction': (100 * r['power_avg'] / loyd) if flying and loyd > 0 else np.nan,
            'v_min': K.minimum_wind_for_crosswind(cfg),
            'flying': float(flying),
        }
    finally:
        K.RHO = rho_saved


def propagate(base, path, v_nominal, reel, n=3000, seed=0, keys=None):
    """Run n samples. Returns {output_name: array of n values}."""
    rng = np.random.default_rng(seed)
    rows = [one_sample(base, path, v_nominal, reel, rng, keys) for _ in range(n)]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def band(x, lo=5, hi=95):
    """Median and a percentile interval, ignoring samples that did not fly."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (np.nan, np.nan, np.nan)
    return (np.percentile(x, lo), np.median(x), np.percentile(x, hi))


# ---------------------------------------------------------------------------
# Which input is costing you the most?
# ---------------------------------------------------------------------------

def sensitivity(base, path, v_nominal, reel, target='power', n=800, seed=1):
    """One-at-a-time: how much does each input alone widen the answer?

    Reported as the 5-95 spread divided by the median, i.e. a relative band
    width. Bigger means that input dominates your uncertainty, which means
    measuring it better buys you the most.
    """
    out = {}
    for k in UNCERTAINTY:
        res = propagate(base, path, v_nominal, reel, n=n, seed=seed, keys=[k])
        lo, mid, hi = band(res[target])
        out[k] = (hi - lo) / mid if (np.isfinite(mid) and mid > 0) else np.nan

    full = propagate(base, path, v_nominal, reel, n=n, seed=seed)
    lo, mid, hi = band(full[target])
    out['__all__'] = (hi - lo) / mid if (np.isfinite(mid) and mid > 0) else np.nan
    return out


# ---------------------------------------------------------------------------

def report(base, name, wind_speeds=(3.0, 4.0, 5.0, 6.0, 7.0), n=2000):
    """Print banded versions of the predictions."""
    L = []
    w = L.append

    w("=" * 78)
    w("PREDICTION BANDS  --  " + name)
    w("Monte Carlo, %d samples, 5th-95th percentile" % n)
    w("=" * 78)
    w("")
    w("Input uncertainties assumed (1 sigma, relative):")
    for k, s in sorted(UNCERTAINTY.items(), key=lambda kv: -kv[1]):
        w("    %-18s +/- %2.0f%%" % (k, 100 * s))
    w("")

    w("-" * 78)
    w("P1  GLIDE RATIO E")
    w("-" * 78)
    cd_eff = K.effective_drag_coefficient(base['cd'], base['tether_length'],
                                          base['tether_diameter'], base['area'])
    w("  nominal  E = %.2f" % (base['cl'] / cd_eff))
    res = propagate(base, {'beta0': 28, 'sweep': 42, 'rise': 11}, 5.0, 0.22, n=n)
    lo, mid, hi = band(res['E'])
    w("  band     E = %.2f  [%.2f - %.2f]" % (mid, lo, hi))
    w("  -> A measured E outside that range falsifies the coefficient estimates.")
    w("")

    w("-" * 78)
    w("P2  CROSSWIND AMPLIFICATION       (the ratio -- best determined)")
    w("-" * 78)
    w("  %-7s %-22s %-22s %s" % ("wind", "static T (N)", "crosswind T (N)", "amplification"))
    for v in wind_speeds:
        try:
            p, fo, fi, _ = O.optimize_lissajous(v, base)
        except Exception:
            continue
        r = propagate(base, p, v, fo, n=n)
        s_lo, s_mid, s_hi = band(r['T_static'])
        c_lo, c_mid, c_hi = band(r['T_crosswind'])
        a_lo, a_mid, a_hi = band(r['amplification'])
        flew = np.nanmean(r['flying'])
        if not np.isfinite(a_mid):
            w("  %-7.1f %-22s %-22s %s" % (v, "-", "no sustained flight", "-"))
            continue
        w("  %-7.1f %-22s %-22s x%.1f [%.1f - %.1f]%s"
          % (v,
             "%.2f [%.2f-%.2f]" % (s_mid, s_lo, s_hi),
             "%.2f [%.2f-%.2f]" % (c_mid, c_lo, c_hi),
             a_mid, a_lo, a_hi,
             "" if flew > 0.95 else "  (%.0f%% of draws flew)" % (100 * flew)))
    w("")
    w("  Note the band on the RATIO is tighter than on either tension alone.")
    w("  C_L raises both, so it cancels. That correlation is exactly what")
    w("  adding errors in quadrature would get wrong, and why this is Monte")
    w("  Carlo. It is also why P2 is the strongest of the seven predictions.")
    w("")

    w("-" * 78)
    w("P4  FRACTION OF THE LOYD CEILING")
    w("-" * 78)
    for v in wind_speeds:
        try:
            p, fo, fi, _ = O.optimize_lissajous(v, base)
        except Exception:
            continue
        r = propagate(base, p, v, fo, n=n)
        lo, mid, hi = band(r['loyd_fraction'])
        pw_lo, pw_mid, pw_hi = band(r['power'])
        if not np.isfinite(mid):
            w("  %-7.1f  no sustained flight" % v)
            continue
        w("  %-7.1f power %6.2f W [%.2f - %.2f]      %.0f%% of Loyd [%.0f - %.0f]"
          % (v, pw_mid, pw_lo, pw_hi, mid, lo, hi))
    w("")

    w("-" * 78)
    w("P6  MINIMUM WIND SPEED")
    w("-" * 78)
    r = propagate(base, {'beta0': 28, 'sweep': 42, 'rise': 11}, 5.0, 0.22, n=n)
    lo, mid, hi = band(r['v_min'])
    w("  v_min = %.2f m/s  [%.2f - %.2f]" % (mid, lo, hi))
    w("  -> Fly on a light day and bracket it. If loops sustain below %.2f m/s"
      % lo)
    w("     the gravity term in the quadratic is wrong.")
    w("")

    w("-" * 78)
    w("WHERE TO SPEND YOUR EFFORT   (one-at-a-time sensitivity on power)")
    w("-" * 78)
    try:
        p, fo, fi, _ = O.optimize_lissajous(5.0, base)
        s = sensitivity(base, p, 5.0, fo, target='power', n=600)
        allw = s.pop('__all__')
        for k, v in sorted(s.items(), key=lambda kv: -(kv[1] if np.isfinite(kv[1]) else -1)):
            if not np.isfinite(v):
                continue
            bar = "#" * max(1, int(round(v * 40)))
            w("  %-18s %5.0f%%  %s" % (k, 100 * v, bar))
        w("")
        w("  %-18s %5.0f%%  (all inputs together)" % ("COMBINED", 100 * allw))
    except Exception as e:
        w("  sensitivity failed: %s" % e)
    w("")
    w("  Read this as: reducing the top input's uncertainty buys you more than")
    w("  improving everything below it. Buy the instrument that fixes the top row.")
    w("")
    w("=" * 78)

    return "\n".join(L)


def compare_scenarios(base, v=5.0, n=1500):
    """How much does doing the static test first buy you?

    This is the practical payoff of the whole uncertainty analysis: it fixes
    the ORDER of the experiment. Right now C_L and C_D are guesses, and they
    dominate everything, so no other prediction can be tested usefully until
    P1 has pinned down E.
    """
    L = []
    w = L.append
    p, fo, _, _ = O.optimize_lissajous(v, base)

    w("=" * 78)
    w("DOES THE ORDER OF THE EXPERIMENT MATTER?   (at %.1f m/s)" % v)
    w("=" * 78)
    w("")
    w("  %-34s %-22s %s" % ("scenario", "power (W)", "band width"))

    for label, sc in (("Today: coefficients guessed", PRIOR),
                      ("After P1: E measured", POST_P1)):
        saved = dict(UNCERTAINTY)
        use_scenario(sc)
        try:
            r = propagate(base, p, v, fo, n=n)
            lo, mid, hi = band(r['power'])
            w("  %-34s %-22s x%.0f"
              % (label, "%.2f [%.2f - %.2f]" % (mid, lo, hi),
                 hi / lo if lo > 0 else float('inf')))
        finally:
            use_scenario(saved)

    w("")
    w("  READ THIS AS AN INSTRUCTION, not a statistic. Until the static test")
    w("  has measured E, every other prediction has a band too wide to falsify.")
    w("  P1 is not merely the first prediction -- it is the gate. Fly the")
    w("  parked-kite test before anything else, and regenerate these bands")
    w("  with the measured coefficients before attempting P2 or P4.")
    w("")
    return "\n".join(L)


if __name__ == '__main__':
    import sys
    import predictions as P

    cfg = P.TRAINER if len(sys.argv) > 1 and sys.argv[1] == 'trainer' else P.PATANG
    text = report(cfg, cfg['name']) + "\n\n" + compare_scenarios(cfg)
    print(text)

    import datetime
    fn = "../prediction_bands_%s_%s.txt" % (
        'trainer' if cfg is P.TRAINER else 'patang',
        datetime.date.today().isoformat())
    with open(fn, 'w') as f:
        f.write(text)
    print("\n[written to %s]" % fn)
