# What this project is, and why — the full explanation

Written to be read start to finish. No prior knowledge assumed beyond school algebra and
trigonometry. If you can solve a quadratic and you know what sine and cosine are, you have enough.

---

## The one-paragraph version

A kite flying *across* the wind pulls dramatically harder than a kite sitting still — hard enough
that people are building power stations out of them. The equations describing this were written down
in 1980 and have been tested on kites from 20 kW up to megawatt scale. **Nobody has tested them on a
kite small enough to buy at a festival stall**, because there's no money in it. This project builds
those equations into working code, commits to seven specific numerical predictions in advance, and
then flies a 12-gram paper patang with a load cell on the line to find out which predictions survive.

---

# Part 1 — Why a kite can make electricity

## The problem with wind turbines

Look at a wind turbine. Big tower, three long blades. Here's what almost nobody realises:

**Almost all the power comes from the outer third of the blades.**

Power from a moving surface depends on how fast it's moving. The blade tips travel at 100+ mph. The
part near the hub barely moves at all. So the inner two-thirds of that enormous blade is essentially
dead weight — it exists only to hold the useful bit out in the air.

And to hold those tips up there you need a tower, a foundation, a hub, a gearbox, and a crane to
install it all. Hundreds of tonnes of steel and concrete, most of it doing nothing but *positioning*.

So somebody asked the obvious question: what if we just fly the tip, and hold it up with a string?

That's airborne wind energy. Replace the tower with a rope.

## Making your own wind

Here's the trick that makes it work, and it's genuinely surprising the first time you see it.

Fly a kite normally. It sits in the sky and tugs gently. Your patang in a 5 m/s breeze pulls with
about **130 grams** of force. Noticeable, not exciting.

Now swoop it fast across the sky, side to side. Something dramatic happens.

Think about sticking your hand out of a car window. At walking pace you feel nothing. At 60 mph it
nearly rips your arm backwards — and the actual wind outside might be dead calm. **You made your own
wind by moving through the air.**

A kite does exactly this. When it flies sideways fast, the air it feels is the real wind *plus* its
own motion. And there's a feedback loop: more speed makes more lift, which makes it go faster still.

## What stops it

Drag. The kite accelerates until drag balances the driving force. Where that balance lands depends
on one number:

$$E = \frac{\text{lift}}{\text{drag}} = \frac{C_L}{C_D}$$

This is the **glide ratio**. It's the same number that tells a glider pilot how far they travel per
metre of height lost — a glide ratio of 40 means 40 km forward per km down. It's a measure of how
good the wing is.

Working through the force balance gives a beautifully simple result:

$$\text{kite speed} \approx E \times \text{wind speed}$$

Your patang has $E \approx 2.3$. In a 5 m/s wind it flies at about 12 m/s — faster than the wind
itself. A good foil kite with $E = 5$ would do 25 m/s. A Makani M600 with $E = 11$ does over 60 m/s.

## Why this matters so much

Aerodynamic force goes as **speed squared**. So if flying crosswind makes the kite go $E$ times
faster, the force goes up by $E^2$.

| Kite | $E$ | $E^2$ |
|---|---|---|
| Patang | 2.3 | 5.4 |
| Foil kite | 4.5 | 20 |
| Makani M600 | 11 | 121 |

**A rigid wing pulls over a hundred times harder flying crosswind than sitting still.** That factor
is the entire reason the airborne wind energy industry exists.

(In practice you don't get the full $E^2$, because the kite isn't pointed straight downwind and
because you're reeling line out. For a real flight it works out around $E^2/2$. For your patang the
model predicts about **×2.8** — smaller than a foil kite's, because $E$ is squared and yours is low.)

## Miles Loyd, 1980

An engineer called Miles Loyd wrote all this down in a six-page paper called *Crosswind Kite Power*.
His result is the ceiling every airborne wind project is measured against:

$$P = \frac{2}{27}\,\rho\, A\, v^3\, C_L\, E^2$$

Decoded:

| Symbol | Meaning | Typical |
|---|---|---|
| $\rho$ | air density | 1.225 kg/m³ |
| $A$ | kite area | 0.12 m² for your patang |
| $v$ | wind speed | 5 m/s |
| $C_L$ | lift coefficient — how "lifty" | ~0.7 |
| $E$ | glide ratio | ~2.3 |

Two things dominate. **Wind cubed** — double the wind, eight times the power. And **glide ratio
squared** — that's the crosswind trick.

### Where the 2/27 comes from

This bit is worth understanding because it's clean calculus and judges like it.

You get power by letting line out against the pull: $P = F \times v_r$, where $v_r$ is your reel-out
speed. But reeling out reduces the wind the kite feels — if you reel out at speed $v_r$ into a wind
$v$, the kite only experiences $(v - v_r)$.

Force goes as the square of that, so:

$$P \propto (v - v_r)^2 \times v_r$$

Write $f = v_r/v$ (reel-out as a fraction of wind speed):

$$P \propto (1-f)^2 f$$

Differentiate and set to zero:

$$\frac{d}{df}\left[(1-f)^2 f\right] = (1-f)^2 - 2f(1-f) = (1-f)(1 - 3f) = 0$$

So $f = 1/3$. **Reel out at one third of the wind speed.** The maximum value is
$(2/3)^2 \times (1/3) = 4/27$, and combined with the $\tfrac12$ from $\tfrac12\rho v^2$ you get
$2/27$.

That's it. The most quoted number in airborne wind energy comes from a one-line derivative.

## The pumping cycle

So you reel line out and it spins a generator. Eventually you run out of line. Then what?

You **depower** the kite — change its angle so it stops pulling hard — and wind the line back in
cheaply. Then power up and go again.

Earn a lot on the way out, pay back a little on the way in. Net profit. It's a yo-yo, and it's
called a **pumping cycle**. SkySails sells commercial systems that work exactly this way.

$$\bar{P} = \frac{P_{\text{out}} t_{\text{out}} - P_{\text{in}} t_{\text{in}}}{t_{\text{out}} + t_{\text{in}}}$$

## Why figure-eights?

Two reasons.

**Practical:** flying in circles winds your tether up like a telephone cord. A figure-eight unwinds
itself — it turns left as much as it turns right.

**Physical:** the kite needs to keep moving crosswind to keep its speed bonus. A stationary kite gets
nothing. So it has to keep sweeping.

*But hold that thought — we found something interesting about this. See Part 7.*

---

# Part 2 — What we are actually asking, and the honest story of how we got there

## What we tried first

The original plan was to find something genuinely new. We generated a hypothesis, then searched the
literature to check it was open. Three times:

**Attempt 1 — "gravity distorts the optimal figure-eight at small scale."** The reasoning: gravity
matters more for small kites, so the optimal loop should be lopsided.

*Dead.* A 2026 paper had already studied it on a 150 m² kite and found exactly that effect. Worse,
my scaling argument was wrong: the gravity number is
$G = mg/(\tfrac12\rho v^2 A C_L)$, and if mass scales with area then **the area cancels** — $G$
doesn't care how big the kite is. Large *rigid* wings are actually the gravity-dominated case.

**Attempt 2 — "Reynolds-number effects change the optimal path at small scale."**

*Dead.* "You can't match Reynolds number in a scaled model" is textbook knowledge in wind energy.
People have been compensating for it for decades.

**Attempt 3 — "tether sag is unmeasured at small scale."** Three separate papers flagged sag as
unresolved in their own work, and one said a small testbed was the right tool for it. Promising!

*Dead.* There's a paper literally titled *The Influence of Tether Sag on Airborne Wind Energy
Generation*. Plus analytical models going back to 2011 and a 2025 paper that includes sag in state
estimation *with experimental data*. And the killer line: **"the straight-tether assumption is
adequate for estimating power generation in a small-scale airborne wind energy system."** We'd have
been measuring something already shown to be negligible exactly where we'd be measuring it.

## What that taught us

Airborne wind energy has had 45 years since Loyd, several well-funded university groups, and serious
commercial R&D. **The open problems left need megawatt hardware, autopilot teams, or CFD clusters.**
None are reachable with a kite and a load cell.

That's not a statement about your ability. It's a statement about what's left in a mature field.

## What we're doing instead, and why it's still real science

Here's the thing those three searches *did* turn up. Look at what professional small-scale campaigns
actually managed to collect, in their own words:

> "The available data covered only **five separate cycles**, which is by far not sufficient for a
> meaningful statistical analysis" — WES 2019, on a 25 m² kite

> Flight B lasted **53 seconds**. Flight A's load cell data was lost entirely to sensor failure.
> — WES 2025, small-scale testbed

Five cycles. Fifty-three seconds. That's the published state of the art at the small end — not
because it's hard, but because their hardware is expensive and their field time is scarce.

**You can beat that in one afternoon, and keep beating it every weekend for a season.**

So the question becomes:

> **How far do published airborne wind energy models hold when you shrink the kite by five orders of
> magnitude, and where exactly do they break?**

That is honest replication at an untested scale. It's not new physics. It *is* real science — nobody
has checked, the models claim to be general, and a general claim that has never been tested at one
end of its range is worth testing.

### Why it's fair to test their models with your kite

You can't compare your *numbers* to a 150 m² kite's numbers — different scale, meaningless
comparison. A judge would catch that immediately.

But those papers publish **equations**, and every one of them claims to be general. Loyd's formula
has no size restriction in it. The tether drag formula has no size restriction. They were derived
from physics, not fitted to big kites.

So: **feed your measurements into their equations and see whether they predict what you measured.**
If a formula claims universality and fails at 0.12 m², that's a finding — and it's a fair test
precisely *because* the formula claimed to be universal.

---

# Part 3 — The model, piece by piece

Everything below is in [`model/kite.py`](model/kite.py).

## The kite lives on a sphere

The tether has a fixed length $L$. So the kite can't go anywhere it likes — it's stuck on the surface
of an invisible sphere of radius $L$, centred on your anchor.

That means it has only **two degrees of freedom** instead of three, and we can describe its position
with two angles:

- $\beta$ (beta) — **elevation**, how high up from the horizon
- $\phi$ (phi) — **azimuth**, how far to the side of straight-downwind

## The single most important angle

Let $\theta$ be the angle between the tether and the wind direction. Some geometry gives:

$$\cos\theta = \cos\beta\,\cos\phi$$

And power scales as $\cos^3\theta$.

**This is why kites fly low and centred.** Straight downwind at low elevation, $\theta = 0$ and
$\cos^3\theta = 1$ — full power. At 60° elevation, $\cos^3(60°) = 0.125$ — you've lost 87% of your
power. The interactive page shades the wind window by $\cos^3\theta$ so you can see it directly.

## How fast does the kite actually fly?

This is the heart of the model. Balance forces *along* the direction of travel: the kite accelerates
until the forward pull from tilting its lift is cancelled by drag, plus whatever component of gravity
is fighting it.

Writing that out and rearranging gives a **quadratic in the flight speed** $v_a$:

$$v_a^2 - \underbrace{E\left(v_w\cos\theta - v_r\right)}_{b}\,v_a + \underbrace{\frac{m g_\tau}{\tfrac12\rho A C_D}}_{c} = 0$$

where $g_\tau$ is how much of gravity opposes the motion (positive when climbing).

Solve it the usual way: $v_a = \frac{b + \sqrt{b^2 - 4c}}{2}$. We take the larger root — the fast,
flying solution.

### And here's the beautiful part

**What if $b^2 - 4c < 0$?**

Then the discriminant is negative and **there is no real solution**. Not a small answer — *no answer
at all*.

Physically: the kite cannot go fast enough to generate the lift it needs to climb against its own
weight. It stalls and falls out of the sky.

Every kite flyer has watched this happen on a light day. **An algebraic condition on a discriminant
predicts something you can see happen in a field.** That's prediction P6, and it's the sharpest thing
in this project.

For your patang the model says this happens below **2.72 m/s**.

## Tether drag

The tether isn't free. It's a long thin cylinder being dragged sideways through the air, and it has
its own drag.

Its speed varies from zero at the anchor to full kite speed at the top. Integrating along its length
and referring the result to the kite's area gives the standard formula everyone uses:

$$C_{D,\text{total}} = C_{D,\text{kite}} + \frac{C_{D,t}\, d\, L}{4A}$$

The factor of 4 comes from that integration.

**For long tethers this is the dominant loss.** For your patang on 50 m of line it's about 17% of
total drag; at 100 m it's 29%.

### And it creates an optimal tether length

Longer tether → more drag (bad), but also higher altitude → stronger wind (good). Since power goes as
$v^3$ and wind grows with height as $v \propto z^\alpha$, there's a genuine trade-off with an
interior maximum. Doing the calculus:

$$L_{\text{opt}} = \frac{3\alpha\, C_{D,\text{kite}}}{\kappa(2 - 3\alpha)}, \qquad \kappa = \frac{C_{D,t}d}{4A}$$

For your patang that's **42–103 m** depending on how rough the ground is. Your 50 m line is already
close to optimal, and you can *test this* by flying at 20 m, 50 m and 100 m and looking for a peak in
the middle. That's P5.

## The turning radius floor

A kite can't turn arbitrarily tightly — it's a wing, not a point. The literature sets the limit at
about **five wingspans**.

This isn't a detail. Without it the optimiser found a cheat: shrink the figure-eight to a *point*
while still collecting the crosswind speed bonus. That describes a kite that isn't actually flying
crosswind, and the quasi-steady formula doesn't notice. The constraint is what makes the model
behave, and it's *active* at the optimum — meaning it genuinely shapes the answer.

## What is deliberately NOT in the model

Being explicit about this matters, because a judge will ask and because it tells you what to expect:

- **Turbulence** — real wind is gusty and gusts don't average out when power goes as $v^3$
- **Unsteady aerodynamics** — the kite has inertia; it can't change speed instantly
- **Canopy deformation** — a paper kite bends and flutters; the model treats it as rigid
- **Control losses** — a real pilot isn't a perfect optimiser

All four of these *cost* power. So **the model should over-predict**, and the size of that gap is the
interesting quantity. We call the lumped correction $\zeta_{\text{turn}}$ and we'll fit it from your
data. Fitting one coefficient from measurement is honest modelling; the fitted value is itself a
result.

---

# Part 4 — The seven predictions, and why each one matters

These live in [`predictions_patang_2026-08-06.txt`](predictions_patang_2026-08-06.txt), dated and not
to be edited after flying.

## Why pre-registration matters

If you measure first and then tune the model to match, you've proven nothing — a flexible enough
model fits anything. If you **commit to numbers in advance**, the model can be *wrong*.

Being capable of being wrong is what makes something science rather than curve-fitting. Judges
recognise this and almost no school project does it.

---

### P1 — Glide ratio from a parked kite ⭐ the gate

Park the kite at the top of the wind window, directly downwind, and hold it still. Three forces act:
lift (up), drag (downwind), weight (down). The tether takes the resultant.

Resolve horizontally and vertically:

$$D = T\cos\beta, \qquad L = T\sin\beta + W$$

$$\boxed{E = \frac{T\sin\beta + W}{T\cos\beta}}$$

**Look at what this needs: a load cell, a phone inclinometer, and a kitchen scale.** No anemometer.
No wind tunnel. No assumed coefficients. No kite area. You measure the single most important number
in the whole project with about ₹400 of equipment.

*Prediction:* $E = 2.32 \pm 0.3$, the same at every wind speed.
*Falsified if:* $E$ drifts systematically with wind speed — which would mean Reynolds effects or
canopy deformation matter at this scale, itself an interesting result.

---

### P2 — The crosswind tension jump ⭐ the headline

Once you know $E$ from P1, the model predicts the tension jump with **zero free parameters**:

$$\frac{T_{\text{crosswind}}}{T_{\text{static}}} \approx E^2\cos^2\theta \times (\text{reel-out term})$$

*Prediction:* **×2.5 to ×3.3** for your patang.

This is the demonstration that sells the project. Your kite pulls 130 g parked, then you start flying
eights and it pulls 400 g. On a bigger kite it's more dramatic still — a 3 m² foil goes from 2.6 kg
to 23 kg. A judge can watch the number jump on your logger.

---

### P3 — The shape of the tension trace

Not just the average — the *waveform*. Within one figure-eight:

$$T(t) \propto \cos^2\theta(t)$$

So tension **peaks at the centre** of the eight (kite lowest and most directly downwind) and **dips
at the two outer edges**.

*Prediction:* **two peaks per full figure-eight**, because the kite crosses the powered centre twice
per loop.
*Falsified if:* one peak per loop, or peaks at the edges.

This is a much stronger test than a single average number, because it tests the *geometry* of the
model, not just its magnitude.

---

### P4 — What fraction of Loyd you achieve

Published full-scale systems reach **10–30%** of Loyd's ceiling.

*Prediction:* your patang achieves **45%**.

**This is the one I'd bet gets falsified, and that makes it the most valuable.** It sits above the
published band. Either the patang genuinely does better (plausible — it's feather-light, so gravity
barely touches it), or the model is missing losses that only show up in reality.

If measurement comes back at 20%, you've found a specific place a published model over-predicts at
small scale — and the gap between 45% and 20% *is* your turning-loss coefficient, measured rather
than assumed.

---

### P5 — An optimal tether length exists

*Prediction:* power per unit area peaks at 42–103 m, not at the longest line you own.
*Falsified if:* power rises monotonically with length.

Test by flying at three lengths on the same day.

---

### P6 — The minimum wind speed ⭐ the sharpest

From the discriminant, as derived above.

*Prediction:* sustained figure-eights become impossible below **2.72 m/s**.
*Falsified if:* loops sustain well below that, or fail well above.

Fly on a light day and bracket it. This is the most falsifiable claim in the project, and it comes
straight out of "does this quadratic have a real root."

---

### P7 — The tether drag formula

*Prediction:* $E$ falls with tether length exactly as $C_L / (C_D + C_{D,t}dL/4A)$ says it should —
from 2.58 at 20 m down to 1.98 at 100 m.
*Falsified if:* $E$ is flat with length (formula overstates tether drag) or falls faster (understates
it).

Just repeat P1 at several line lengths.

---

# Part 5 — The optimiser, and why those specific algorithms

In [`model/optimize.py`](model/optimize.py). The algorithm choices aren't arbitrary — they follow
from one fact about the problem.

## The objective function has a cliff in it

Remember P6: where the discriminant goes negative, the kite stalls and power is **undefined**. Not
zero, not small — undefined.

So the function we're maximising has a vertical wall in it. That single fact rules out most standard
optimisers:

| Method | Verdict |
|---|---|
| Gradient-based (SLSQP, BFGS) | ❌ computes garbage finite-difference gradients across the cliff |
| Nelder–Mead, Powell | ✅ derivative-free — only ever *compares* values, doesn't care |
| Differential evolution | ✅ derivative-free *and* global, survives a multi-peaked landscape |

**Being able to explain why you rejected gradient methods is worth more to a judge than the result
itself.** It shows you chose rather than copied.

## Three stages

1. **Coarse grid sweep — always first.** Yes it's the dumb method. Run it anyway: it shows the whole
   shape of the objective, catches multiple peaks before a local method fools you, and makes a good
   poster figure.
2. **Nelder–Mead** to polish the best grid point.
3. **Differential evolution** for the 15-dimensional version, where grid search is hopeless (10
   values per axis would be $10^{15}$ evaluations).

## Two tricks worth knowing

**Penalty instead of undefined.** Rather than returning "no answer" on stall, we return a number that
gets *worse* the deeper into stall you are. That turns a vertical wall into a slope the optimiser can
walk down. It's called an exterior penalty method.

**Reparameterisation instead of constraints.** Instead of checking "is elevation between 15° and 70°"
we *change variables* so it's impossible to be outside that range. Constraints that can't be violated
never need checking.

## The correctness test that caught a real bug

We describe the flight path two ways:

- **Lissajous** — 3 parameters, always a figure-eight
- **Fourier series** — 15 coefficients, any closed loop at all

And here's the key structural fact: **a Lissajous figure-eight is a special case of the Fourier
path** (set $a_1 = $ sweep, $c_2 = $ rise, everything else zero).

So the Fourier optimum **must** be at least as good as the Lissajous optimum. It's searching a bigger
space that contains the smaller one. If it ever does worse, the optimiser is broken.

**On the first run, it did worse.** By 40%.

The cause was real: my Fourier code divided the elevation coefficients by the number of harmonics,
capping the "rise" parameter at 4° when the optimum needs 12°. **The Lissajous path wasn't reachable
inside my supposedly more general search space.** Fixed, round-trip tested, and the global search is
now warm-started from the Lissajous solution so the property holds by construction.

That test cost about ten lines and caught a bug that would have produced a wrong headline result.
Build tests that *can* fail.

## What it found

Give the optimiser fifteen free coefficients and no constraint on shape, and it returns:

```
a = [44.3, 0.1, -0.2]     c = [0.0, 10.8, 0.1]     d ≈ 0
```

Pure first harmonic in azimuth, pure second harmonic in elevation — **that is exactly a Lissajous
figure-eight**, with everything else set to zero.

**The classical figure-eight is within 5% of optimal**, and that 5% comes from retuning elevation and
reel-out speed, not from a better shape. That's a real answer to "why does everyone fly
figure-eights," and you *found* it rather than assuming it.

---

# Part 6 — Uncertainty, and why it changed the experiment

In [`model/uncertainty.py`](model/uncertainty.py).

## Why a single number is a bad prediction

"The power will be 2.32 W" **cannot be wrong.** Any measurement differs from it, so "close enough"
becomes a judgement call and the comparison proves nothing.

"The power will be between 1.0 and 4.9 W" **can be wrong.** If you measure 8 W, the model is
falsified.

So we draw every uncertain input from a distribution a few thousand times, push each draw through the
model, and report the 5th–95th percentile. That's Monte Carlo.

## The error I made, which the code caught

I assumed P2 (the tension *ratio*) would be better determined than either tension alone — a high
$C_L$ raises both, so surely it cancels?

It doesn't:

$$\text{static} \sim C_L, \qquad \text{crosswind} \sim C_L E^2 = \frac{C_L^3}{C_D^2}, \qquad \text{ratio} \sim \frac{C_L^2}{C_D^2}$$

Nothing cancels. The ratio depends on the **square** of an already-uncertain glide ratio, making it
*more* sensitive to the coefficients, not less.

The real argument for Monte Carlo isn't that it beats adding errors in quadrature. It's that **it
stops you reasoning your way to a confident wrong answer.**

## The finding that changed the plan

Vary one input at a time and see how much the answer moves:

```
cl                231%  ████████████████████████████████████
cd                210%  ████████████████████████████████
wind_speed        193%  █████████████████████████████
tether_diameter    26%  ████
area               22%  ███
rho                10%  ██
mass                0%
```

The top three dominate by an order of magnitude — and two of them are currently **guesses**. Which
gives:

| | Predicted power at 5 m/s | Band width |
|---|---|---|
| Today, coefficients estimated | 2.30 W `[0.17 – 11.89]` | **×68** |
| After P1 measures $E$ | 2.33 W `[1.04 – 4.85]` | **×5** |

A ×68 band can't falsify anything.

**So P1 isn't just first on the list — it's a gate.** Do the parked-kite test, put the measured $E$
into the code, regenerate, and only *then* attempt P2 or P4. The uncertainty analysis told us the
order of the experiment, which is not what I expected it to do.

Note also that wind speed is third and it's the one you *can't* fix by measuring the kite. A handheld
anemometer at head height while the kite is 30 m up is exactly the ±20% assumed here. Only a logged
mast fixes that.

---

# Part 7 — The design chart, and the thing it uncovered

In [`model/design_chart.py`](model/design_chart.py), figures in [`figures/`](figures/).

The chart sweeps kite area against wind speed, solving a full trajectory optimisation at every point,
and answers: *given this kite and this wind, what do you get?* It spans your 0.12 m² patang to
Makani's 54 m² M600 — a factor of 450 in area and 250,000 in power.

## The bug that wasn't a bug

The first version showed **nothing above 1 m² able to fly.** Obviously wrong — the M600 flew fine.

Except the model was right, and telling me my configuration was impossible:

- The five-wingspan turning limit grows as **√area**
- I'd scaled tether length as $A^{0.16}$ — far too weakly

So big kites were being asked to fly loops tighter than physically possible. **Tether length has to
scale with wingspan.** Fixed.

## And then the interesting part

While debugging I measured how tight a turn each path shape can hold:

| Path | Tightest turn it can manage |
|---|---|
| Lissajous figure-eight | **0.18 × tether length** |
| Circle | **0.22 × tether length** |

Why? A figure-eight **reverses direction** at each azimuth extreme, and that reversal is a
near-cusp — a very tight turn. A circle never reverses; its curvature is uniform.

**So past a certain size the figure-eight is geometrically impossible and only the circle survives.**

That is very likely why Makani flew circles rather than eights. It's also exactly the "circular
versus figure-of-eight" debate in the 2026 paper — and our model reproduces the reason independently.

*(One trap, in case you present this: a circle must have its azimuth amplitude divided by
$\cos\beta_0$. Without that you trace a circle in **angle space**, which on a sphere is an ellipse
with a tight end — defeating the whole point. My first attempt scored worse than the figure-eight
because of exactly this.)*

## The honest caveat

The model says circles beat figure-eights **everywhere** on pure power. That's not a contradiction of
practice — it's a statement about what's in the model. **There's no penalty here for winding the
tether up**, and tether twist is precisely why the field flies eights.

The advantage of the figure-eight is operational, not aerodynamic. Say that plainly and it's a
strength; hide it and it's a hole.

## The non-dimensional version

The second figure removes size entirely. Everything collapses onto two numbers:

- $E$ — the glide ratio
- $G = mg/(\tfrac12\rho v^2 A C_L)$ — the **gravity number**, how heavy the kite is relative to the
  lift it can make

Any two kites sharing $(E, G)$ land on the same point *regardless of size*. On the chart, the 1.5 m²
trainer and the 3 m² foil sit on exactly the same spot despite one being twice the area.

**This is what licenses a 0.12 m² patang to say anything about a 54 m² Makani wing**, and it's the
single most important idea for defending the project's scope.

Two readings:

- **$E$ dominates.** Contours run almost vertically. Power coefficient tracks glide ratio and barely
  notices gravity until $G > 0.5$. Another reason P1 is the gate.
- **Your patang survives light wind because it's flimsy.** At 12 grams its $G$ is 0.06 — an order of
  magnitude below a power foil's 0.43. The "serious" kites are the ones that fall out of the sky in
  light wind.

---

# Part 8 — What you actually do in a field

## Equipment

| Item | Spec | Approx. |
|---|---|---|
| Load cell | **1 kg** straight bar — static pull is only 39–262 g | ₹250 |
| Amplifier | HX711 24-bit | ₹120 |
| Logger | ESP32 or Arduino Nano + microSD | ₹400 |
| Line angle | phone inclinometer app, or a BNO085 IMU at the anchor | ₹0–700 |
| Anemometer | cup type on a pole, logged | later |

Note the load cell is **1 kg, not 200 kg** — your patang pulls under half a kilogram even flying
crosswind. That's a ₹250 part.

## The sequence

1. **Measure the kite.** Area of the diamond, mass on a kitchen scale, line length. Put the real
   numbers into `PATANG` in `model/predictions.py`, regenerate, and commit — *before* flying, so the
   pre-registration is honest.
2. **P1 first, always.** Park the kite, log tension and angle at several wind speeds. Compute $E$.
3. **Update the model** with the measured $E$ and regenerate the bands. They should tighten by ~13×.
4. **P7 while you're there** — repeat P1 at 20 m, 50 m, 100 m of line.
5. **Then P2 and P3.** Fly figure-eights, log the tension trace, look for the jump and the two peaks
   per loop.
6. **P6 on a light day.** Find the wind speed where sustained loops become impossible.
7. **P5** across three tether lengths on one day.

Fly often. Remember what you're competing with: five cycles and 53 seconds.

## Safety

Your patang pulls under half a kilogram — it genuinely cannot hurt you, which is one good reason to
start there. But:

- **Never use manja.** Glass-coated line maims birds and people every year and is banned in several
  states. Plain cotton or polyester only. It also gives you a known line diameter, which the tether
  drag term needs.
- Open ground, nobody downwind.
- If you later move to a trainer foil pulling 9 kg, that changes — quick-release at the anchor,
  gloves, never wrap line around your hands.
- Check local airspace rules, and **never fly near a BARC facility.**

---

# Part 9 — What's still to do

- [ ] **You:** measure the patang and regenerate the predictions
- [ ] **You:** verify your GitHub email so Pages can build (the runner is being refused)
- [ ] Load cell wiring guide
- [ ] Mine Makani's 1000-page *Energy Kite Report* for stated limitations
- [ ] Extend the model to fly-gen (Loyd's "drag mode") so Makani's public flight logs become
      comparable — their data is free on BigQuery and covers a 600 kW machine

---

# Glossary

| Term | Meaning |
|---|---|
| **AWE** | Airborne wind energy |
| **Glide ratio $E$** | Lift ÷ drag. How good the wing is. Power goes as $E^2$ |
| **$C_L$, $C_D$** | Lift and drag coefficients — dimensionless measures of lift and drag |
| **Crosswind flight** | Flying across the wind to make your own headwind |
| **Loyd limit** | The 1980 theoretical ceiling, $\frac{2}{27}\rho A v^3 C_L E^2$ |
| **Pumping cycle** | Reel out under load, depower, reel in cheap, repeat |
| **Reel-out factor $f$** | Reel speed ÷ wind speed. Optimal is 1/3 |
| **Elevation $\beta$** | Angle up from the horizon |
| **Azimuth $\phi$** | Angle sideways from straight-downwind |
| **$\theta$** | Angle between tether and wind. Power goes as $\cos^3\theta$ |
| **Gravity number $G$** | $mg/(\tfrac12\rho v^2 A C_L)$ — weight relative to lift |
| **Quasi-steady** | Assumes forces balance at every instant; ignores inertia |
| **Lissajous curve** | A figure-eight made by combining two oscillations at a 1:2 frequency ratio |
| **Discriminant** | The $b^2 - 4ac$ in a quadratic. Negative means no real solution |
| **Pre-registration** | Writing predictions down before measuring, so they can be wrong |
| **Patang** | Indian festival fighter kite |
