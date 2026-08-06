"""
Mathematical model of a tethered kite flying crosswind figure-eights.

Everything here is "quasi-steady": we assume that at every instant the kite
has settled into force balance. That is an approximation, but it is the
standard one in airborne-wind-energy literature and it keeps the maths
readable.

Coordinate system
-----------------
Ground anchor at the origin. Wind blows along +x. z is up.

The kite is stuck on the surface of a sphere of radius L (the tether length),
so it only has two degrees of freedom:

    beta  = elevation angle, measured up from the horizon
    phi   = azimuth angle, measured sideways from straight-downwind

Unit vector from anchor to kite:

    e_r = (cos(beta)cos(phi),  cos(beta)sin(phi),  sin(beta))

The single most important geometric quantity is theta, the angle between
the tether and the wind:

    cos(theta) = e_r . x_hat = cos(beta) * cos(phi)

Power scales as cos^3(theta), so this one number does a lot of work.
"""

import numpy as np

RHO = 1.225   # air density, kg/m^3 (sea level, 15 C)
G = 9.81      # gravity, m/s^2


# ---------------------------------------------------------------------------
# 1. Tether drag  ->  effective glide ratio
# ---------------------------------------------------------------------------

def effective_drag_coefficient(cd_kite, tether_length, tether_diameter,
                               kite_area, cd_tether=1.0):
    """Total drag coefficient of kite + tether, referred to the kite's area.

    The tether is a long cylinder being dragged sideways through the air.
    Its speed varies linearly from 0 at the anchor to the full kite speed at
    the top. Integrating v^2 along its length and dividing by the kite's
    reference speed gives the classic factor of 1/4:

        extra C_D = C_D,tether * d * L / (4 * A)

    This is usually the DOMINANT loss for long tethers, and it is what
    creates an optimal tether length (see optimal_tether_length).
    """
    return cd_kite + cd_tether * tether_diameter * tether_length / (4.0 * kite_area)


def effective_glide_ratio(cl, cd_kite, tether_length, tether_diameter,
                          kite_area, cd_tether=1.0):
    """E_eff = C_L / C_D,total. Power goes as E^2, so this matters enormously."""
    cd = effective_drag_coefficient(cd_kite, tether_length, tether_diameter,
                                    kite_area, cd_tether)
    return cl / cd


# ---------------------------------------------------------------------------
# 2. Static kite  -- this is what the cheap validation experiment measures
# ---------------------------------------------------------------------------

def glide_ratio_from_static_test(tension, beta_deg, kite_mass):
    """Measure E = L/D from a PARKED kite using only a load cell and a phone.

    Park the kite at the top of the wind window, directly downwind, and let it
    sit still. Three forces act on it: lift (up, perpendicular to the wind),
    drag (downwind), and weight (down). The tether takes the resultant.

    Resolving horizontally and vertically:

        D = T * cos(beta)
        L = T * sin(beta) + W

    so
        E = L/D = (T sin(beta) + W) / (T cos(beta))

    Note what is NOT needed: no anemometer, no kite area, no air density,
    no assumed C_L. This is the highest-value measurement in the project,
    because every power prediction scales as E^2.
    """
    b = np.radians(beta_deg)
    weight = kite_mass * G
    return (tension * np.sin(b) + weight) / (tension * np.cos(b))


def static_tension_and_angle(wind_speed, kite_area, cl, cd, kite_mass):
    """Forward prediction: where does a parked kite sit, and how hard does it pull?

    Inverse of glide_ratio_from_static_test. Used to check the model against
    the first, easiest field measurement.
    """
    q = 0.5 * RHO * wind_speed ** 2 * kite_area   # dynamic pressure * area
    lift = q * cl
    drag = q * cd
    weight = kite_mass * G
    beta = np.arctan2(lift - weight, drag)
    tension = np.hypot(lift - weight, drag)
    return tension, np.degrees(beta)


# ---------------------------------------------------------------------------
# 3. Crosswind flight speed  --  the quadratic, including gravity
# ---------------------------------------------------------------------------

def apparent_wind_speed(wind_speed, cos_theta, reel_speed, e_eff, cd_eff,
                        kite_area, kite_mass, climb_component):
    """Solve for how fast the kite flies through the air.

    Force balance ALONG the flight path. The kite accelerates until the
    forward pull it gets from tilting its lift vector is cancelled by drag
    plus the component of gravity opposing its motion:

        (1/2) rho A C_L v_a (v_w cos(theta) - v_r)  =  (1/2) rho A C_D v_a^2  +  m g_tau

    Rearranged, this is a quadratic in v_a:

        v_a^2  -  b v_a  +  c  =  0

        b = E_eff * (v_w cos(theta) - v_r)      "driving term"
        c = m g_tau / ((1/2) rho A C_D_eff)     "gravity penalty"

    We take the larger root (the fast, flying solution).

    KEY RESULT: if b^2 < 4c the discriminant is negative and there is NO real
    solution. Physically, the kite cannot sustain that climb -- it stalls and
    falls out of the sky. This predicts a minimum wind speed for crosswind
    flight, which is a real thing every kite flyer has seen in light wind, and
    it is much more severe for small kites (large m/A) than for big ones.

    Parameters
    ----------
    climb_component : float
        g_tau / g, i.e. the vertical component of the kite's direction of
        travel. Positive when climbing, negative when diving.

    Returns
    -------
    float or nan
        Apparent wind speed, or nan if the kite cannot fly here.
    """
    b = e_eff * (wind_speed * cos_theta - reel_speed)
    if b <= 0:
        return np.nan                      # kite is being outrun by the winch

    c = kite_mass * G * climb_component / (0.5 * RHO * kite_area * cd_eff)

    disc = b * b - 4.0 * c
    if disc < 0:
        return np.nan                      # gravity wins -- kite stalls
    return 0.5 * (b + np.sqrt(disc))


def tether_tension(v_a, kite_area, cl, kite_mass, beta):
    """Radial force balance: aerodynamic lift minus the radial bit of weight."""
    lift = 0.5 * RHO * kite_area * cl * v_a ** 2
    return max(lift - kite_mass * G * np.sin(beta), 0.0)


# ---------------------------------------------------------------------------
# 4. The figure-eight path  --  a Lissajous curve
# ---------------------------------------------------------------------------

def lissajous_state(s, params):
    """Three-parameter figure-eight.

        phi(s)  = sweep * sin(s)
        beta(s) = beta0 + rise * sin(2s)

    Because beta oscillates at exactly twice the frequency of phi, the path
    closes into a figure-eight that crosses itself at the centre. Frequency
    ratio 1:2 is what makes it an eight rather than an oval; this is the
    standard Lissajous result.

    Real kites fly figure-eights rather than circles for a completely
    practical reason: a circle winds the tether up, an eight does not.

    Returns (beta, phi, dbeta/ds, dphi/ds), all in radians.
    """
    sweep = np.radians(params['sweep'])
    rise = np.radians(params['rise'])
    beta0 = np.radians(params['beta0'])

    phi = sweep * np.sin(s)
    beta = beta0 + rise * np.sin(2.0 * s)
    dphi = sweep * np.cos(s)
    dbeta = 2.0 * rise * np.cos(2.0 * s)
    return beta, phi, dbeta, dphi


def fourier_state(s, params):
    """General closed periodic path as a truncated Fourier series.

        phi(s)  =        sum_k  a_k sin(ks) + b_k cos(ks)
        beta(s) = beta0 + sum_k  c_k sin(ks) + d_k cos(ks)

    The Lissajous figure-eight is the special case a_1 = sweep, c_2 = rise,
    everything else zero. That nesting is deliberate: the Fourier optimum
    MUST beat the Lissajous optimum, because the Lissajous path lies inside
    its search space. If it ever does worse, the optimizer is broken -- a
    free correctness test.

    Coefficients are in degrees, harmonics k = 1..N.
    """
    beta0 = np.radians(params['beta0'])
    a = np.radians(np.asarray(params['a']))     # phi, sine terms
    b = np.radians(np.asarray(params['b']))     # phi, cosine terms
    c = np.radians(np.asarray(params['c']))     # beta, sine terms
    d = np.radians(np.asarray(params['d']))     # beta, cosine terms

    k = np.arange(1, len(a) + 1)
    ks = np.outer(k, np.atleast_1d(s))          # (n_harmonics, n_samples)

    phi = (a[:, None] * np.sin(ks) + b[:, None] * np.cos(ks)).sum(axis=0)
    beta = beta0 + (c[:, None] * np.sin(ks) + d[:, None] * np.cos(ks)).sum(axis=0)
    dphi = (k[:, None] * (a[:, None] * np.cos(ks)
                          - b[:, None] * np.sin(ks))).sum(axis=0)
    dbeta = (k[:, None] * (c[:, None] * np.cos(ks)
                           - d[:, None] * np.sin(ks))).sum(axis=0)
    return beta, phi, dbeta, dphi


def _spectral_derivative(f, order=1):
    """Exact derivative of a periodic function sampled on a uniform grid.

    Differentiate in Fourier space: multiply the k-th mode by (i k)^order.
    For a BAND-LIMITED periodic function this is exact, not approximate --
    and our paths are band-limited by construction, because they ARE finite
    Fourier series. Same reason the trapezoid rule is the right integrator
    here: periodicity makes the spectral method the natural one.
    """
    n = len(f)
    k = np.fft.fftfreq(n, d=1.0 / n)          # integer wavenumbers
    return np.real(np.fft.ifft(np.fft.fft(f) * (1j * k) ** order))


def path_geometry(s, params, state_fn=lissajous_state):
    """Vectorized path geometry over a full loop.

    `s` must be a uniform grid over [0, 2*pi). Returns, as arrays:

        cos_theta     angle between tether and wind (the cos^3 term lives here)
        beta          elevation
        speed_factor  |d e_r / ds|, arc length per unit s
        climb         vertical component of the direction of travel
        radius        LOCAL RADIUS OF CURVATURE of the flight path, in metres

    The radius matters because a kite cannot turn arbitrarily tightly. Without
    that constraint the optimizer discovers a degenerate solution: shrink the
    loop to a point, keep the crosswind speed bonus, and collect infinite
    power from a kite that is not actually flying crosswind. Eijkelhof et al.
    (WES 2024) set the lower bound at about 5 wingspans.
    """
    beta, phi, _, _ = state_fn(s, params)

    e = np.stack([np.cos(beta) * np.cos(phi),
                  np.cos(beta) * np.sin(phi),
                  np.sin(beta)])                       # (3, n) unit sphere

    de = np.stack([_spectral_derivative(e[i], 1) for i in range(3)])
    dde = np.stack([_spectral_derivative(e[i], 2) for i in range(3)])

    speed_factor = np.linalg.norm(de, axis=0)
    safe = np.maximum(speed_factor, 1e-12)

    climb = de[2] / safe
    cos_theta = np.cos(beta) * np.cos(phi)

    # Curvature of r(s) = L * e(s):  kappa = |r' x r''| / |r'|^3
    cross = np.cross(de.T, dde.T).T
    kappa = np.linalg.norm(cross, axis=0) / safe ** 3
    radius = np.where(kappa > 1e-12, 1.0 / np.maximum(kappa, 1e-12), 1e9)

    return cos_theta, beta, speed_factor, climb, radius


# ---------------------------------------------------------------------------
# 5. Cycle-averaged power
# ---------------------------------------------------------------------------

def reel_out_power(wind_speed, kite, path, reel_factor, n=256,
                   state_fn=lissajous_state):
    """Average power during the reel-out (power-generating) phase.

    Integrates around one full figure-eight. Two things worth noting:

    1. We average over TIME, not over the path parameter s. The kite spends
       different amounts of time in different parts of the loop and that
       matters a lot.

    2. The integration is a plain uniform-spacing trapezoid sum, which looks
       naive but is not: for a SMOOTH PERIODIC integrand over a full period,
       the trapezoid rule converges EXPONENTIALLY rather than at the usual
       O(h^2). Simpson's rule is actually worse here. Our integrand is
       periodic in s by construction, so the simplest method is the best one.
       (Trefethen & Weideman, SIAM Review 2014.)

    Returns (mean_power_W, mean_tension_N, loop_period_s, stalled_fraction).
    """
    L = kite['tether_length']
    e_eff = effective_glide_ratio(kite['cl'], kite['cd'], L,
                                  kite['tether_diameter'], kite['area'])
    cd_eff = effective_drag_coefficient(kite['cd'], L,
                                        kite['tether_diameter'], kite['area'])
    v_r = reel_factor * wind_speed
    zeta = kite.get('zeta_turn', 1.0)   # lumped turning loss, FITTED from data

    s_vals = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    ds = 2.0 * np.pi / n

    cos_theta, beta, speed_factor, climb, radius_rel = path_geometry(
        s_vals, path, state_fn)
    turn_radius = radius_rel * L        # path_geometry works on the unit sphere

    # --- flight-speed quadratic:  v_a^2 - b v_a + c = 0 -------------------
    b = e_eff * (wind_speed * cos_theta - v_r)
    c = kite['mass'] * G * climb / (0.5 * RHO * kite['area'] * cd_eff)
    disc = b * b - 4.0 * c

    flying = (b > 0) & (disc >= 0)      # real root exists => kite can fly here
    v_a = np.where(flying, 0.5 * (b + np.sqrt(np.maximum(disc, 0.0))), np.nan)

    lift = 0.5 * RHO * kite['area'] * kite['cl'] * np.nan_to_num(v_a) ** 2
    T = np.maximum(lift - kite['mass'] * G * np.sin(beta), 0.0)

    dt = np.where(flying, L * speed_factor * ds / np.where(flying, v_a, 1.0), 0.0)

    total_time = dt.sum()
    if total_time <= 0:
        return 0.0, 0.0, np.inf, 1.0, 0.0

    return (zeta * (T * v_r * dt).sum() / total_time,
            (T * dt).sum() / total_time,
            total_time,
            1.0 - flying.mean(),
            float(np.min(turn_radius)))


def reel_in_power(wind_speed, kite, beta_in_deg, reel_in_factor):
    """Cost of hauling the kite back in.

    Trick of the trade: fly the kite to high elevation (small cos(theta)) and
    DEPOWER it (drop C_L by changing angle of attack) before winching in. Both
    cut the force you have to pull against.

    Reeling in makes the apparent wind LARGER, not smaller, because you are
    dragging the kite into the wind:  v_a = v_w cos(theta) + v_in.
    """
    b = np.radians(beta_in_deg)
    cos_theta = np.cos(b)                      # phi = 0, parked straight downwind
    v_in = reel_in_factor * wind_speed

    v_a = wind_speed * cos_theta + v_in
    force = 0.5 * RHO * kite['area'] * np.hypot(kite['cl_depowered'],
                                                kite['cd']) * v_a ** 2
    return force * v_in, force


def cycle_power(wind_speed, kite, path, reel_factor, reel_in_factor,
                beta_in_deg, stroke_length, state_fn=lissajous_state):
    """Full pumping cycle: earn on the way out, pay on the way back.

        P_avg = (P_out * t_out  -  P_in * t_in) / (t_out + t_in)

    This is the number the whole project is trying to maximize.
    """
    p_out, t_mean, _, stalled, r_turn = reel_out_power(
        wind_speed, kite, path, reel_factor, state_fn=state_fn)
    p_in, f_in = reel_in_power(wind_speed, kite, beta_in_deg, reel_in_factor)

    t_out = stroke_length / (reel_factor * wind_speed)
    t_in = stroke_length / (reel_in_factor * wind_speed)

    p_avg = (p_out * t_out - p_in * t_in) / (t_out + t_in)
    duty = t_out / (t_out + t_in)

    return {
        'power_avg': p_avg,
        'power_out': p_out,
        'power_in': p_in,
        'tension_mean': t_mean,
        'force_reel_in': f_in,
        't_out': t_out,
        't_in': t_in,
        'duty_cycle': duty,
        'stalled_fraction': stalled,
        'min_turn_radius': r_turn,
    }


# ---------------------------------------------------------------------------
# 6. Benchmarks and closed-form results
# ---------------------------------------------------------------------------

def loyd_limit(wind_speed, kite_area, cl, e_ratio):
    """Loyd 1980, the theoretical ceiling.

        P = (2/27) rho A v^3 C_L (C_L/C_D)^2

    Derivation of the 2/27: tension goes as (v_w - v_r)^2 and power is
    tension * v_r, so P ~ (1-f)^2 f with f = v_r/v_w. Differentiating,
    d/df[(1-f)^2 f] = 0 gives f = 1/3 and a maximum value of 4/27. Combined
    with the leading 1/2 from (1/2) rho v^2, that is 2/27.

    Real systems achieve 10-30% of this. The gap is what the rest of this
    model is about.
    """
    return (2.0 / 27.0) * RHO * kite_area * wind_speed ** 3 * cl * e_ratio ** 2


def optimal_tether_length(cd_kite, tether_diameter, kite_area,
                          shear_exponent=0.14, cd_tether=1.0):
    """A longer tether reaches stronger wind but adds drag. There is an optimum.

    Wind grows with height as v ~ z^alpha, and altitude z = L sin(beta), so

        P  ~  L^(3 alpha) / (C_D,kite + kappa L)^2 ,     kappa = C_D,t d / (4A)

    Setting dP/dL = 0 gives

        L_opt = 3 alpha C_D,kite / (kappa (2 - 3 alpha))

    A pleasant little calculus result with a directly testable prediction:
    fly at three tether lengths and the power should peak in the middle.
    """
    kappa = cd_tether * tether_diameter / (4.0 * kite_area)
    a = shear_exponent
    return 3.0 * a * cd_kite / (kappa * (2.0 - 3.0 * a))


def wind_at_height(v_ref, z_ref, z, roughness=0.03):
    """Log wind profile. Wind is much faster up high, and power goes as v^3,
    so this correction is larger than every other one in the model combined.

        v(z) = v_ref * ln(z/z0) / ln(z_ref/z0)

    roughness z0: 0.01-0.03 m open grass, 0.05 m crops, 0.3-1.0 m suburban.
    """
    return v_ref * np.log(z / roughness) / np.log(z_ref / roughness)


def minimum_wind_for_crosswind(kite, cos_theta=0.85, max_climb=0.5):
    """Smallest wind speed at which the kite can fly the whole figure-eight.

    Solves b^2 = 4c for v_w, i.e. the point where the discriminant of the
    flight-speed quadratic hits zero at the steepest climbing point of the loop.
    With reel-out at the ideal f = 1/3:

        E_eff * v_w * (cos(theta) - 1/3) = 2 sqrt(m g climb / ((1/2) rho A C_D))
    """
    e_eff = effective_glide_ratio(kite['cl'], kite['cd'], kite['tether_length'],
                                  kite['tether_diameter'], kite['area'])
    cd_eff = effective_drag_coefficient(kite['cd'], kite['tether_length'],
                                        kite['tether_diameter'], kite['area'])
    c = kite['mass'] * G * max_climb / (0.5 * RHO * kite['area'] * cd_eff)
    return 2.0 * np.sqrt(c) / (e_eff * (cos_theta - 1.0 / 3.0))
