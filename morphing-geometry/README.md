# Morphing wing shape generator (geometry team)

Takes a baseline airfoil plus a few morph settings and writes morphed airfoil
coordinate files. Every shape is also checked automatically, to see whether the
structure could actually build it.

- Baselines:
  - **NASA SC(2)-0714** supercritical airfoil (**default**), the transonic
    case, loaded from NASA's published coordinates (see
    [below](#nasa-sc2-0714-transonic-case)).
  - **NACA 64(1)-412**, the original baseline from the handoff, built from its
    two parts: the 64(1)-012 thickness form and the a = 1.0 mean line with
    design lift coefficient 0.4. Run any script with `config/naca641412.yaml`.
  - Any other airfoil given as a `.dat` coordinate file.
- Morph: the shape changes only ahead of the front spar and behind the
  hinge, which sits at the rear spar or anywhere behind it. The wing box
  between the spars stays fixed. There are two methods:
  - **v1 `camber`:** shift the camber line, keep the thickness. Simple, but the
    skins have to stretch.
  - **v2 `constant_length`:** the skins bend but keep their length.
- Output: `.dat` files for XFOIL and CFD meshing, XYZ text files in mm for
  SolidWorks, PNG plots, and a `summary.csv` for sweeps.
- Studies: the **deflection range** (how far the TE can droop and reflex
  before a check fails) and the **rear spar position** (how far back the spar
  can go and still reach the range you need). See
  [below](#hinge-deflection-range-and-rear-spar-position).

Written in Python 3.12 with numpy, scipy, matplotlib and pyyaml; tests use pytest.

---

## Quick start

**Setup (once).** On Krishna's Mac the `morphgeom` conda env already exists, in
the native ARM64 miniforge at `~/miniforge3`. Don't use the old x86
`~/opt/miniconda3`. On another machine, create the env from inside the
project folder:

```bash
conda env create -f environment.yml
```

**Each new terminal**, inside the project folder:

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate morphgeom
```

Then run any script:

| Command | What it does |
|---|---|
| `python scripts/make_baseline.py` | Unmorphed SC(2)-0714 → `outputs/sc20714/baseline/`. Prints how far it is from the published coordinates. |
| `python scripts/make_shape.py` | One morphed shape → `outputs/sc20714/shapes/` (files, 3 plots, check report). Edit the settings at the top, including `cfg.morph.method`. |
| `python scripts/sweep.py` | Every combination in the `sweep:` settings (both methods by default) → `outputs/sc20714/sweep/` plus `summary.csv`. |
| `python scripts/make_step.py <name>_cad_mm.txt --span <mm>` | Turns any CAD point file into STEP files that SolidWorks opens directly: a flat 2D face, plus a 3D solid if you give a span. |
| `python scripts/deflection_range.py` | Largest droop and reflex that pass every check, for the settings at the top of the script → `outputs/sc20714/range/` (the two limit shapes, a CSV and a plot). About 5 seconds. |
| `python scripts/spar_study.py` | Deflection range at a series of rear spar stations, and the furthest-aft station that still reaches the needed range → `outputs/sc20714/spar_study/`. About 35 seconds. |
| `python scripts/variable_study.py` | Changes one setting at a time (hinge, shape exponent, skin, chord, strain limit, ...) and finds the range each time, then ranks which settings matter most → `outputs/sc20714/variable_study/`. A few minutes. |
| `pytest` | Runs the 27 tests in `tests/` (about 45 seconds). |

With no case file, every script works on the **SC(2)-0714**. Each script also
takes a **case file** for another airfoil or setup, for example the NACA
64(1)-412 (its files go to `outputs/naca641412/`):

```bash
python scripts/make_baseline.py config/naca641412.yaml
```

All settings live in **`config/default.yaml`**. Anything marked
`PLACEHOLDER` there is not a design decision. It is only there so the code runs.
A case file (like `config/naca641412.yaml`) lists only what differs from
`default.yaml`. A misspelt setting name in a case file stops with an error
instead of being ignored.

From your own Python script or notebook:

```python
from morphgeom import load_config, make_shape, write_outputs
cfg = load_config()
cfg.morph.method = "constant_length"
cfg.morph.delta_LE = 0.01
res = make_shape(cfg)            # res.checks has every check value, res.row the summary
write_outputs(res, cfg, "outputs/my_case")
```

(Scripts outside `scripts/` need the project folder on `sys.path`, or run them from the project folder.)

---

## Conventions

- Chord = 1. x runs from 0 (leading edge) to 1 (trailing edge), y is up.
- Points: cosine spacing (bunched at both edges), 200 per surface by default. Leading and trailing edge are included.
- Point order: **Selig** (TE → upper surface → LE → lower surface → TE). Set `points: order: reverse` for the other direction.
- Sign of the morph: `delta > 0` is droop (edge down, more camber). `delta < 0` is reflex (edge up).
- **Morphed shapes are not re-normalised.** The wing box stays where it is and the drooped TE simply sits lower (y = −0.020 for `delta_TE = 0.02`). For CFD, keep the angle of attack relative to the baseline x-axis so every shape is compared at the same wing attitude.
- Summary numbers are printed to 3 decimals, in units that make 3 decimals useful (mm, %, degrees, 1/m). Coordinate files keep 6 decimals.

### File names

`naca641412_TE+0.020_xh0.650_n3.dat`: TE drop 0.020 c, rear spar at 0.650 c, exponent 3.
An LE part (`LE+0.010_xle0.150`) appears only when `delta_LE` is not 0.
A rear spar ahead of the hinge adds `_rs0.600` after the hinge (`TE+0.020_xh0.700_rs0.600_n3`).
v2 shapes end in `_constlen`.
No morph at all gives `naca641412_baseline` (the same for both methods).

### CAD file (SolidWorks)

`<name>_cad_mm.txt` holds three tab-separated columns, X Y Z in **mm**, with Z = 0 and no header.
It is the same single loop as the `.dat`, scaled by `cfg.structure.chord_m`.
To use it: set the part units to mm, then Insert → Curve → Curve Through XYZ Points → Browse.
For a closed trailing edge the first and last points are the same, so the curve closes.
CATIA and most other CAD tools can read the same three-column text.

### STEP files (open directly in SolidWorks)

Instead of rebuilding the profile from points in SolidWorks, you can let the
script make the CAD file:

```bash
python scripts/make_step.py outputs/sc20714_tunnel/baseline/sc20714_baseline_cad_mm.txt --span 203.2
```

Next to the point file it writes:
- `<name>_2d.step`: the closed profile as a flat face at z = 0 (for 2D CFD meshing or sketches).
- `<name>_3d_span<S>mm.step`: a solid, the same face extruded S mm and centred on z = 0. It is only written if you give `--span`.

In SolidWorks: **File → Open**, set the file type to **STEP (*.step; *.stp)**,
and pick the file. It opens as a new part, in mm, with no other steps.

- The outline is one smooth curve through every point, the same as Curve
  Through XYZ Points. A blunt trailing edge is closed with a straight line.
- After writing, the script reads each file back and prints its size, its
  area and its volume. The area matches the area inside the points to within
  0.01%.
- The span is not something this project decides, so there is no default.
  Pass it each time.
- It also works on several files at once (for example `outputs/sc20714/shapes/*_cad_mm.txt`).
- This needs the `pythonocc-core` package (the OpenCascade CAD kernel). It is
  in `environment.yml` and already installed in the `morphgeom` env. Nothing
  else in the project needs it.

---

## How a shape is built

**1. Baseline.** The half-thickness `yt(x)` comes from the 64(1)-012 table in `data/`.
The camber line is the a = 1.0 formula:

```
yc/c   = -(cl_i / (4*pi)) * [ (1 - x) ln(1 - x) + x ln(x) ]
dyc/dx =  (cl_i / (4*pi)) * ln( (1 - x) / x )
```

The thickness is laid **perpendicular to the camber line**, the standard NACA way
(`build` in `morphgeom/surfaces.py`).

**2a. v1 morph, `camber`** (`camber_delta` in `morphgeom/morph.py`).
- Behind the hinge `x_h`: `s = (x − x_h)/(1 − x_h)`, `dy = −delta_TE · s^n`.
- Ahead of the front spar `x_LE`: `s = (x_LE − x)/x_LE`, `dy = −delta_LE · s^n`.
- From the front spar to the hinge, `dy` is exactly zero. By default the
  hinge *is* the rear spar (`x_rear_spar: null`). Set `x_rear_spar` ahead of
  `x_h` to leave fixed skin between the spar and the hinge.
- The new camber line goes through the same construction with the same thickness.
- n = 2: the camber line's slope matches at the spar. n = 3 (default): its curvature matches too.

The report also gives each morph as a plain flap angle: `atan(delta_TE / (1 − x_h))`.

**2b. v2 morph, `constant_length`** (`constant_length` in `morphgeom/morph.py`).
Each skin is cut into many short straight pieces. The morph only *turns*
pieces, never stretches them, so a skin's length can't change. How much a
piece turns grows smoothly from the spar, and `n` keeps its meaning (2 = slope
matches at the spar, 3 = curvature matches too).

- **Trailing edge: "sliding TE".**
  - You can't droop a closed, sharp TE whose skins can't stretch while ribs
    hold them apart: the outer skin needs more length on a bend than the
    inner one. Real morphing TEs let something give. Here it's the lower skin.
  - The **upper skin** runs unbroken to the tip and keeps its length. It bends
    with the same curvature law as v1, just enough to drop the TE by
    `delta_TE`. The TE then moves forward a few mm.
  - The **ribs** (the line from each upper skin point to the lower skin point
    at the same station) keep their length and their angle to the upper skin.
    So the thickness doesn't change and the skins can't cross.
  - The **lower skin** slides a few mm into the TE. The checks report how far.
- **Leading edge: constant-length droop nose** (as in the Swansea and DLR
  morphing-LE work).
  - The skin wrapping round the nose is clamped at both front-spar points and
    keeps its total length.
  - v1's droop shape can't be used as it is, because it needs about 1% more
    skin for a 1% chord droop. So v2 makes **a new droop shape with the same
    skin length**.
  - The nose drops exactly as far as v1's, measured as the mid-height between
    the surfaces at 1% chord. Within that, the shape stays as close to the v1
    droop as it can while bending smoothly.
  - The result: the skin bends mostly just ahead of the front spar, and the
    nose swings down nearly as a whole.
  - Settings:
    - `cfg.morph.le_smoothing` (default 0.01) trades closeness to v1 against
      smoothness.
    - `cfg.morph.le_full_droop = False` drops the droop condition. You then get
      the closest shape to v1, which reaches only about half the droop.

---

## The checks (every shape, value + PASS/FAIL)

The checks use their own dense points: 2001 per surface, with a point exactly at
each spar. So results don't depend on how many points go into the output files.

| # | Check | What it measures | Limit |
|---|---|---|---|
| 1 | Wing box | Largest movement of any skin point inside the box, front spar to rear spar (mm). Points *by station* must be identical; the code checks they are, bit for bit. | `tol_wingbox_mm` (placeholder) |
| 2 | Skin length | % length change of each skin piece: TE upper, TE lower, LE upper, LE lower. The whole LE wrap (spar → nose → spar) and both surface totals are reported too. Also given in mm. In v2 the TE lower skin is designed to slide, so it is reported as the slide, not judged. | `tol_arc_pct` (placeholder) |
| 3 | Bending strain | ε = (t/2)·Δκ. Δκ is the curvature change between matching skin points, matched by arc-length fraction within each piece. Needs `chord_m` and `t_skin_mm`. | `tol_strain_pct` (placeholder) |
| 4 | Smooth joins | Jump in slope (deg) and curvature (1/m) at the front spar and the hinge, on the curve the morph is built on: the camber line (v1), or the upper TE skin and LE skin (v2). Jumps on the other skins are reported as info. | numerical zero |
| 5 | Sanity | Thickness ≥ 0; the loop never crosses itself; TE gap unchanged; no new inflection points (wiggles) in the LE skin. | — |
| info | | Nose radius before/after, nose drop at 1% chord, TE movement, largest thickness change, rough lift change (`dcl_estimate`, see [below](#rough-lift-estimate)). v2 also gives the lower-skin slide and the distance from the v1 shape. | |

Worked example (test 5): a 2.000 mm skin bent from flat to a 0.500 m radius gives
ε = 0.001 m × 2.000 1/m = 0.002 (0.200%).

"Numerical zero" for the curvature jump is 0.01 1/m. A step that small changes a
2 mm skin's strain by less than 0.001%, which shows as 0.000 at 3 decimals. A
real n = 2 step is 0.3–8 1/m.

---

## What the runs show

These are for the NACA 64(1)-412 (`config/naca641412.yaml`, files in
`outputs/naca641412/`), with the placeholder spars (0.15, 0.65), a 1 m chord
and a 2 mm skin. The SC(2)-0714 results are in [its own section](#nasa-sc2-0714-transonic-case).

**Baseline matches published coordinates.**
- With the Abbott & von Doenhoff table, the largest difference from the PDAS
  64(1)-412 table is **0.006% chord** (0.059 mm on 1 m). That is the size of
  the difference between the two thickness tables themselves.
- With PDAS's own thickness table, the difference is 0.000% chord: only their
  rounding is left. That confirms the mean line and the construction.

**v1 against v2**, n = 3:

| | v1 TE 2% | v2 TE 2% | v1 LE 1% | v2 LE 1% |
|---|---|---|---|---|
| Skin length change (worst piece) | +1.284% (TE upper) | 0.000% | +3.774% (LE upper) | 0.000% |
| Lower TE skin slides into TE | – | 3.7 mm | – | – |
| Peak bending strain | 0.078% | 0.100% | 1.013% | 0.308% |
| Nose drop at 1% chord | – | – | 6.8 mm | 6.8 mm |
| Distance from the v1 shape | – | 1.5 mm | – | 5.0 mm |
| Largest thickness change | 0 | 0 | 0 | 0.95 mm |
| TE moves forward | 0 | 4.7 mm | – | – |
| All checks | FAIL (skin length) | **PASS** | FAIL | FAIL (LE hollow only) |

What this means:
- **The v1 camber morph is not length-preserving.** A 2% TE droop stretches
  the upper TE skin by 1.3%. A 1% LE droop stretches the upper nose skin by
  3.8% and squeezes the lower by 2%.
- **v2 fixes the TE cleanly.** The shape stays within 1.5 mm of v1 and every
  check passes. The cost is a lower skin that slides 3.7 mm into the TE, and a
  TE that moves 4.7 mm forward.
- **v2 reaches the full LE droop with the same skin length**, and a third of
  v1's strain. It gets there with a different shape, up to 5 mm from v1's:
  - the skin bends mostly just ahead of the front spar (a gentle shoulder on
    top, a gentle hollow underneath),
  - the nose itself keeps its shape (radius 9.8 mm before and after) and
    swings down.
  - CFD should judge whether that shape keeps the aerodynamic benefit of the
    droop.
- **The only remaining LE flag is the lower-surface hollow** (2 new inflection
  points). v1 has one too at 1% droop. In v2 it is where the skin bends.
  Whether a gentle concave patch ahead of the spar is acceptable is a call for
  the CFD and structures teams.
- **Across other cases** (droop −1% to +2%, front spar at 15% or 25% chord, n = 2
  or 3), v2 always reaches v1's nose drop exactly and nothing crosses. With
  n = 3 the strain is 3–4× lower than v1. A front spar at 25% chord makes it
  easier: a 1% droop needs only 0.13% strain.
- **The follower skins need one order more smoothness.** In v1 the skins follow
  the camber line at a distance, and in v2 the lower TE skin follows the ribs.
  So their curvature only matches exactly at the spar with n = 4. With n = 3
  the leftover step is 0.013–0.027 1/m, about 0.003% strain. It's shown as
  info.

---

## Hinge, deflection range and rear spar position

### The hinge

`x_h` is the **hinge**: where the TE starts to bend. Everything behind it is
free to morph. `x_rear_spar` is the back of the wing box.
- `x_rear_spar: null` (default): the hinge is at the rear spar, so the TE
  bends right from the spar. This gives the morph the most skin to bend.
- `x_rear_spar` ahead of `x_h` (e.g. spar 0.60, hinge 0.70): the skin between
  them stays fixed, for example to make room for actuators. The wing box
  check covers the spars; the smoothness check moves to the hinge.
- A spar behind the hinge stops with an error.

Set them in `make_shape.py`, `deflection_range.py`, a case file, or the
`sweep:` lists (`x_rear_spar: [null, 0.60]`).

### Deflection range: `scripts/deflection_range.py`

For the settings at the top of the script (method, hinge, n), it finds the
**largest droop and the largest reflex that still pass every check**.
- Deflections are given as the **equivalent flap angle**: the angle of a
  straight flap hinged at `x_h` with the same tip drop. This lets you compare
  them with ordinary flaps. With `n = 3` the bending is concentrated towards
  the tip, so right at the TE the surface is about three times as steep as
  this angle suggests.
- **How it searches:** it tries 1°, 2°, 3°, … until a check fails, then halves
  the gap between the last pass and the first fail down to 0.01°. Each answer
  is an angle that *passed*. The search settings are under `range:` in
  `config/default.yaml`.
- `edge = "LE"` does the same for the nose (equivalent angle
  `atan(delta_LE / x_LE)`).
- **Output** (`outputs/<case>/range/`):
  - the two limit shapes as `.dat` and CAD files;
  - a CSV of every shape tried, with every check value;
  - a plot with the limit shapes, the limiting check against angle, and the
    rough lift change.

Results with the placeholders (1 m chord, 2 mm skin, 0.5% strain limit,
hinge at 0.65, n = 3, v2):

| | Largest droop | Largest reflex | Stopped by |
|---|---|---|---|
| SC(2)-0714 | +13.898° (tip down 86.6 mm) | −16.984° (tip up 106.9 mm) | bending strain |
| NACA 64(1)-412 | +14.578° (tip down 91.0 mm) | −16.180° (tip up 101.6 mm) | bending strain |

- **The range depends on the placeholders.** Strain = (skin thickness / 2) ×
  curvature change. So doubling the allowed strain, or halving the skin
  thickness, roughly doubles the range. A longer chord gives the same skin
  more room to bend.
- **v1 (`camber`) stops at about 0.2–0.3°.** Its skins stretch, and they pass
  the 0.1% skin-length limit almost at once. That is why v2 exists. To see
  where v1 would stop if stretching were allowed, set
  `cfg.range.ignore = ["skin_length"]`.
- **The nose (v2) stops at about ±1°,** at the LE shape check: the gentle
  hollow under the nose (open question 11).
- These are **structural** limits. Whether the flow stays attached at the
  ends of the range is for XFOIL or CFD (next steps).

### Rear spar position: `scripts/spar_study.py`

**The trade-off.** A rear spar further back makes a bigger wing box, but
leaves less skin behind it to bend. So the morphing TE reaches a smaller
angle before the skin is over-strained, and each degree gives less lift.

**What the script does.** It runs the deflection-range search at each station
in `spar_study.x_rear_spar`, with the hinge at the spar. Then it reports the
furthest-aft station that still reaches `need_droop_deg` and
`need_reflex_deg` (PLACEHOLDERS: +10° and −5°, open question 14). For each
station it also gives:
- the wing box length;
- the rear spar height (skin to skin);
- the box's share of the section area;
- the rough lift per degree.

SC(2)-0714 with the placeholders (front spar 0.15, v2, n = 3):

| Rear spar x/c | Spar height (mm) | Box share of area | Largest droop | Largest reflex | Lift per degree (est.) | Reaches +10° / −5° |
|---|---|---|---|---|---|---|
| 0.50 | 132.568 | 50.089% | 17.891° | −22.586° | 0.158 | yes |
| 0.60 | 113.878 | 63.353% | 15.352° | −18.984° | 0.143 | yes |
| 0.65 | 100.012 | 69.067% | 13.898° | −16.984° | 0.135 | yes |
| 0.70 | 84.085 | 73.984% | 12.328° | −14.867° | 0.126 | yes |
| 0.75 | 66.700 | 78.010% | 10.633° | −12.633° | 0.116 | yes |
| 0.80 | 48.877 | 81.090% | 8.820° | −10.305° | 0.104 | no |
| 0.85 | 32.707 | 83.254% | 6.891° | −7.906° | 0.091 | no |

- The **droop need decides it** here. The furthest-aft station that reaches it
  is 0.75; between the tested stations, about **x = 0.767**.
- For the 64(1)-412 the same placeholders give about 0.781.
- Treat this as the *method*, not the answer. The number moves with the
  needed range, the strain limit, the skin thickness and the chord, which
  are all placeholders.
- The CSV and plot are in `outputs/<case>/spar_study/`. For a finer look, set
  a closer list of stations at the top of the script.

### Which settings matter: `scripts/variable_study.py`

**What it does.** It changes **one setting at a time**, with everything else
left at its default, and runs the deflection-range search at each value.
Then it ranks the settings by how much they move the range. The settings and
values are listed under `variable_study:` in `config/default.yaml`; add,
remove or change lines there. The values are **test values**, chosen wide
enough to show which settings matter. They are not design choices.

```bash
python scripts/variable_study.py
```

It takes a few minutes. It writes to `outputs/sc20714/variable_study/`:
- `sc20714_TE_variable_study_constlen.csv`, and `_LE_` for the leading edge:
  one row per value tried;
- the matching `.png`: the ranking at the top, then one small plot per
  setting, all on the same angle scale.

**TE results, SC(2)-0714** (v2, placeholders otherwise; needed range +10° / −5°):

| Setting | Values tried | Largest droop over them | Largest reflex over them | Effect |
|---|---|---|---|---|
| strain limit | 0.25 – 1.5 % | +7.7° to +22.2° | −8.6° to −30.0° | biggest |
| chord | 0.5 – 2 m | +7.7° to +20.9° | −8.6° to −29.0° | biggest |
| skin thickness | 4 – 1 mm | +7.7° to +20.9° | −8.6° to −29.0° | biggest |
| shape exponent n | 2 – 5 | +0.1° to +16.0° | −0.1° to −19.8° | big: n = 2 fails, best at 2.8 – 3 |
| hinge x/c | 0.80 – 0.55 | +8.8° to +16.7° | −10.3° to −20.9° | big |
| solid TE strip | 0 – 3% chord | +13.6° to +14.4° | −16.5° to −17.9° | small |
| front spar x/c | 0.10 – 0.25 | +13.9° | −17.0° | none |

What this shows:

1. **Chord, skin thickness and strain limit act as ONE number.** The strain
   check only sees

   *skin thickness ÷ (2 × strain limit) ÷ chord*

   which is the tightest extra bend the skin may take, measured in chords.
   It is 0.200 at the defaults (2 mm, 0.5%, 1 m). Runs with the same value
   give exactly the same range, whichever of the three was changed (the
   last plot panel shows them on one curve).
   - So the structures team only has to agree this one number.
   - The +10° droop need is met when it is about 0.3 or less. For example,
     3 mm of skin on a 1 m chord at 0.5% just misses, at +9.98°.
2. **Hinge.** Each 0.05 c further aft costs about 1.5° of droop and 2° of
   reflex, as in the spar study.
   - The TE range depends on the hinge only. Where the rear spar sits ahead
     of the hinge makes no difference (a test checks this).
3. **Shape exponent n.**
   - Higher n (4, 5) concentrates the bend towards the tip, which roughly
     halves the range.
   - n = 2 has a real curvature jump at the hinge and stops at about 0.1°.
   - Between 2 and 3 the joint is smooth, but v2's lower skin gets a strain
     spike right at the hinge. The spike grows the closer you look: for
     n = 2.5 it is 0.28% with 2,001 check points, 0.43% with 8,001 and 0.88%
     with 32,001. The curvature-join check rightly stops n = 2.5 at about 3°.
   - **n = 2.8** passes every check, and its result is the same with 16×
     finer check points: **+16.0° / −19.8°**, about 15% more than n = 3.
     It still has a weak form of that spike, so ask structures before using
     it (open question 16). n = 3 stays the default.
4. **Solid TE strip.** This is the width at the tip left out of the strain
   check.
   - It hardly matters on the SC(2)-0714, whose TE is blunt.
   - On the 64(1)-412, whose TE is sharp, it matters a lot. With no strip,
     the lower skin's corner at the tip shows 11% strain at only 1°, and the
     range collapses to 0.2°. With 0.5% chord it is back to +14.1° / −15.6°.
5. **Front spar.** No effect on the TE range.
6. **Strain limit 1.5%.** This is the one place the search stops for another
   reason (hollow points on the plot). The v2 solver cannot go past about
   22° of droop, and the reflex reaches the 30° search cap.

**LE results, SC(2)-0714.** The LE range is only about ±1°, and it is set by
the "no new wiggles" check (`le_shape`), not by strain. So skin thickness and
LE smoothing make no difference. Some settings roughly double it:
- n = 4 or 5: about ±2°;
- front spar at 0.10: about ±1.9°;
- `le_full_droop: false`: about +1.8° / −1.8°.

So the LE limit is a question about the shape rule itself (open question 11),
not about the skin. `n_modes` = 2 cannot be solved, so the study tries 3 to 6.

**64(1)-412** (`python scripts/variable_study.py config/naca641412.yaml`):
- the ranking is the same, apart from the solid TE strip (point 4);
- n = 2.8 gives +16.7° / −18.7° (n = 3: +14.6° / −16.2°);
- its default range is +14.6° / −16.2°.

Each setting was changed on its own, so the study does not show what happens
when two change at once (for example a thinner skin *and* a hinge further
back). Run that combination as a case file with `deflection_range.py`.

### Rough lift estimate

`aero.py` estimates the lift change of a morph at the same angle of attack,
using thin-airfoil theory. That means a thin airfoil in smooth, low-speed
flow, with the flow following the surface.
- It is checked against two known answers:
  - the 64(1)-412 mean line gives cl = 0.400 at zero angle (its design value);
  - a straight flap matches the textbook flap formula exactly.
- It is fine for **comparing** shapes and spar stations at small deflections.
  That is why the studies report lift **per degree**.
- At large deflections it reads far too high, because the real flow
  separates. Don't quote lift values at the limit angles.
- XFOIL or CFD give the real numbers. Transonic effects (SC(2)-0714) are not
  included at all.

---

## NASA SC(2)-0714 (transonic case)

A 14%-thick supercritical airfoil with design lift coefficient 0.7 and a blunt
trailing edge (gap 0.70% of chord). NASA tested it at Mach 0.6–0.76 in the
Langley 0.3-m Transonic Cryogenic Tunnel (NASA TM-4044), which makes it a good
case for validating transonic CFD. On a swept wing the section sees roughly
M × cos(sweep): an airliner at Mach 0.85 with 25° sweep puts it near 0.77.

**Run it.** It is the default airfoil, so no case file is needed:

```bash
python scripts/make_baseline.py
```

```bash
python scripts/make_shape.py
```

```bash
python scripts/sweep.py
```

Everything goes to `outputs/sc20714/`. (Adding `config/sc20714.yaml` still
works; it gives the same result.)
Coordinates: `data/downloads/sc20714.dat` (NASA TP-2969 design coordinates, via
the UIUC database).

**Wind-tunnel scale, for CFD validation (phase 2).** `config/sc20714_tunnel.yaml`
draws the same airfoil at the NASA test model's chord, **152.4 mm** (6 in).
NASA's Reynolds numbers are based on that chord.

```bash
python scripts/make_baseline.py config/sc20714_tunnel.yaml
```

That writes `outputs/sc20714_tunnel/baseline/sc20714_baseline_cad_mm.txt`.
The ready-made STEP files are in the same folder (open with File → Open in
SolidWorks):
- `sc20714_baseline_2d.step`: the 2D profile as a flat face, which is most
  likely what the CFD team wants.
- `sc20714_baseline_3d_span203.2mm.step`: a solid, 152.4 mm chord × 203.2 mm
  span, from z = −101.6 to +101.6 mm. The NASA model spanned the 203.2 mm
  (8 in) width of the tunnel's 2D test section.

They were made with:

```bash
python scripts/make_step.py outputs/sc20714_tunnel/baseline/sc20714_baseline_cad_mm.txt --span 203.2
```

To draft it by hand instead: Curve Through XYZ Points (part in mm), a
straight line across the 1.067 mm blunt trailing edge, then Planar Surface
(2D) or Extrude (3D).

Use this case for the unmorphed
baseline only. At this chord the strain check would describe a 6-inch test
model, not the wing, so morphing work stays on `config/sc20714.yaml` with the
full-scale chord.

**How it is loaded** (`split_dat` in `morphgeom/surfaces.py`):
- **NASA's frame is kept.** The file has its nose at (0, 0) and its TE
  midpoint at (1, −0.013). Re-normalising it to put the TE on y = 0 would tilt
  it by 0.74° and shift every angle of attack against the NASA data.
- **The published points are smoothed.** NASA printed these coordinates to 5
  decimals, but they scatter like 4-decimal values (about 0.013 mm per metre
  of chord). A curve forced through every point has curvature that zig-zags
  at every point (about ±1.5 per chord). That would swamp the strain and
  wiggle checks and show as small pressure wiggles in CFD.
  - So camber and thickness are fitted as smoothly as possible while staying
    within that scatter. The scatter is measured from the data itself.
  - The result: rms 0.020 mm and at most 0.050 mm from any published point,
    on a 1 m chord. The curvature is smooth and the LE radius is 3.03% of
    chord. The nose point, TE height and TE gap match the file exactly.
- The baseline check plot (`outputs/sc20714/baseline/sc20714_baseline_check.png`)
  shows the built shape over the file's points.

**First results** (same placeholder spars, chord and skin as above, n = 3):

| | v1 TE 2% | v2 TE 2% | v1 LE 1% | v2 LE 1% |
|---|---|---|---|---|
| Skin length change (worst piece) | +1.625% (TE upper) | 0.000% | +2.784% (LE upper) | 0.000% |
| Lower TE skin slides into TE | – | 4.8 mm | – | – |
| Peak bending strain | 0.067% | 0.101% | 0.909% | 0.256% |
| Nose drop at 1% chord | – | – | 8.1 mm | 8.1 mm |
| Distance from the v1 shape | – | 3.1 mm | – | 4.1 mm |
| TE moves forward | 0 | 6.2 mm | – | – |
| All checks | FAIL (skin length) | **PASS** | FAIL (skin length, strain) | FAIL (LE hollow only) |

- The pattern matches the 64(1)-412. The thicker SC(2) section makes the TE
  skins stretch more under v1 (1.6% against 1.3%), so the v2 lower skin slides
  further (4.8 mm).
- As on the 64(1)-412, v2's LE droop leaves a gentle hollow in the lower
  surface just ahead of the front spar.

**Blunt trailing edge.** Both morphs keep the 7 mm gap, turning it with the TE.
The `.dat` file and the CAD file are then *open* loops: they start at the upper
TE corner and end at the lower one. In SolidWorks, draw a straight line across
the TE to close the profile (`make_step.py` adds that line for you). The CFD team should say how they want the blunt
base meshed (open question 5).

---

## Choices made in the code (and why)

1. **Smooth thickness between table points.** `yt = sqrt(x)·g(x)`, with `g` a cubic
   spline. This keeps the round nose well behaved, and slope and curvature are
   continuous.
   - The Abbott & von Doenhoff table has no LE radius, so the fit extrapolates
     one: 1.055% chord. The same method on the PDAS table gives 1.002% against
     its published 1.007%.
   - If someone looks up the book's LE radius, add `# le_radius = ...` to the
     table file.
2. **Nose (x < 0.5% chord).** The a = 1.0 mean-line slope goes to infinity at
   the LE. Laying thickness along it puts a microscopic kink in the nose. NACA
   tables have no points between 0 and 0.5%. So ahead of 0.5% chord, the
   thickness direction continues smoothly from its value there, and the nose
   curvature falls off steadily from the tip. Published stations (0.5% and up)
   are unaffected.
3. **Trailing-edge strip.** The a = 1.0 camber line's own curvature also grows
   without limit at the TE, in roughly the last 1% of chord. The strain check
   skips `te_exclude` = 1% of chord there. It is a placeholder for the width of
   the solid TE strip.
4. **v2 TE uses a sliding lower skin.** A closed TE with two skins that can't
   stretch can only droop if something else gives. Two other options were
   tried first:
   - Letting both skins keep their length while the section changed shape.
     That gave skins that crossed near the tip, or bumps of 8–12 mm.
   - Letting the thickness change freely. Same problems.

   The sliding-TE layout keeps the thickness exact and the shape smooth.
   *Whether a sliding TE suits the structure is a question for the structures
   team* (open question 10).
5. **v2 LE droop shape: bend the skin, don't reshape the outline.** The
   obvious way to get a same-length droop is to take v1's shape and pull the
   nose back a little (1.7 mm) until the length matches. That was tried. It
   keeps the length, but the skin then has to slide about 5 mm round the
   nose. Material from the tightly curved tip ends up where the skin is much
   flatter, giving about **7.6% strain**: unbuildable. v2 instead bends the
   skin through a few smooth modes, so each piece of skin stays near its own
   place on the nose. That gives 0.3% strain for the same droop.
6. **v2 accuracy.** v2 is built from many short straight pieces, so it carries
   a small discretisation error. It solves on a grid twice as fine as the
   checks use. Stations are placed exactly at the spars without making the
   spacing uneven. The skin lengths in the checks agree to about 1e-6%.
7. **Loaded `.dat` files** (`baseline: source: dat` plus `dat_file:`, for
   example SC(2)-0714 or RAE 2822).
   - **Vertical split.** The split into camber and thickness is approximate.
     It uses a **vertical** split, and these airfoils are rebuilt with
     thickness laid vertically, not perpendicular. A morphed loaded airfoil
     therefore keeps its thickness in the vertical direction. A perpendicular
     split was tried and dropped: near the nose, a tiny camber change rotates
     the cut lines enough that the fit runs away.
   - **Frame.** A file that already has its nose at (0, 0) and its TE at
     x = 1 keeps its own frame. Anything else is shifted, turned and scaled
     to that form, and the report says by how much.
   - **Smoothing.** Camber and thickness are fitted as smoothly as the
     scatter in the published points allows. The scatter is measured from
     the data, and never taken as less than the rounding of the printed
     decimals. The nose point and the TE are pinned to the file.
     `split_dat(..., smooth=False)` passes exactly through every point
     instead.
   - **Fit quality.** `af.split_error` and `af.split_rms` report how far the
     rebuilt shape is from the file's points. Examples: SC(2)-0714 0.050 mm
     at most; the PDAS 64(1)-412 table 0.001 mm.
8. **Curvature wiggles at table points.** The thickness comes from a table
   rounded to 3 decimals at 5% steps. So the curvature plot shows small corners
   at the table stations (most visible at x = 0.05 and 0.40, inside the fixed
   box). They do not affect the morph checks.
   - If CFD shows small pressure wiggles there, generate a dense 64(1)-012 table
     with the PDAS `naca456` program and point `baseline: thickness_file:` at it.
     No code change is needed.

---

## Open questions (answers go into `config/default.yaml`)

1. ~~Python or MATLAB?~~ Python. (A MATLAB version was built first. The Python
   port reproduces its 34 coordinate files and its full `summary.csv` exactly.)
2. ~~SolidWorks or CATIA?~~ SolidWorks (the XYZ file also works for CATIA).
3. Front and rear spar locations, `x_LE` and `x_rear_spar`, and the hinge `x_h`. Placeholders: 0.15 and 0.65, with the hinge at the rear spar. `scripts/spar_study.py` shows how far back the rear spar can go.
4. Skin thickness and chord for the strain check. Placeholders: 2 mm and 1 m. Together with the strain limit (question 6) they only matter through skin thickness ÷ (2 × strain limit) ÷ chord (`variable_study.py`), so agreeing that one number is enough for the range.
5. Coordinate format with the CFD team: point order, points per surface, TE handling. The SC(2)-0714 has a real blunt TE (7 mm on a 1 m chord): how should its base be closed and meshed?
6. Pass/fail limits for skin length, strain and wing-box movement. Placeholders: 0.1%, 0.5%, 0.01 mm.
7. Which flight conditions and lift coefficients the sweep should target. The sweep ranges are placeholders.
8. Width of the solid TE strip left out of the strain check. Placeholder: 1% chord. It barely matters on the blunt SC(2)-0714 TE, but a sharp TE (64(1)-412) needs at least about 0.5% chord (`variable_study.py`).
9. Does the CFD team agree to keep morphed shapes in the baseline frame (not re-normalised)?
10. *New (v2):* is a sliding lower skin at the TE acceptable? If not, which gives: skin stretch (v1) or section thickness?
11. *New (v2):* is the v2 LE droop shape acceptable (it bends just ahead of the front spar, with a gentle hollow underneath)? CFD should compare it with v1's shape at the same droop.
12. ~~*New (SC(2)-0714):* does it replace the 6-series as the main baseline?~~ Yes: it is now the default. The 64(1)-412 stays available through `config/naca641412.yaml`. (The 6-series idea is about tracking the laminar drag bucket; on a supercritical section the gain is shock and aft-loading control.)
13. *New (SC(2)-0714):* do its 14% thickness and design lift coefficient of 0.7 match the wing station you care about? The SC(2) family (NASA TP-2969, all in the UIUC database) has other combinations, e.g. -0410, -0412, -0610, -0612, -0710, -0712. Any of them runs the same way: copy `config/sc20714.yaml` and change the file and name.
14. *New (deflection range):* what range must the morphing TE reach? `spar_study.need_droop_deg` and `need_reflex_deg`. Placeholders: +10° and −5°. Better still as a lift change, once aero/CFD say what the wing needs.
15. *New (hinge):* does the TE start bending right at the rear spar, or is there fixed structure (actuators, a false spar) between them? Set `x_rear_spar` ahead of `x_h` for the second case.
16. *New (shape exponent):* is n a little below 3 acceptable? n = 2.8 passes every check and reaches about 15% more range than n = 3 (+16.0° / −19.8°), but below 3 the v2 lower skin has a (weak, for 2.8) strain spike right at the hinge. Placeholder: n = 3.

## Possible next steps

- CFD of v1 against v2 at the same droop, to see what the constant-length nose
  shape costs or gains aerodynamically.
- Dense thickness table from PDAS `naca456` (choice 8 above).
- CFD of the SC(2)-0714 baseline against NASA TM-4044 (phase 2 gate), using
  the 152.4 mm geometry in `outputs/sc20714_tunnel/baseline/`.
- RAE 2822 through the `.dat` loader, the same way as SC(2)-0714.
- XFOIL or CFD on the two limit shapes from `deflection_range.py`, to check
  the flow stays attached at the ends of the range (the geometry checks
  can't tell).
- Two settings at a time, e.g. skin thickness with hinge position, once the
  structures team narrows the ranges (`variable_study.py` changes one at a time).

---

## Folder layout

```
environment.yml             conda env (python 3.12, numpy, scipy, matplotlib, pyyaml, pytest, pythonocc-core)
config/default.yaml         every setting; placeholders marked. Default airfoil: NASA SC(2)-0714
config/naca641412.yaml      case file for the NACA 64(1)-412 (only what differs from default.yaml)
config/sc20714.yaml         SC(2)-0714 case file: now the same as the default; a template for other .dat airfoils
config/sc20714_tunnel.yaml  SC(2)-0714 at the NASA test model's 152.4 mm chord (CFD validation geometry)
data/                       thickness tables and published reference coordinates (sources inside each file)
  README_reference_files.txt
  downloads/n64012.dat      raw UIUC file the default thickness table came from
  downloads/sc20714.dat     raw UIUC file: NASA SC(2)-0714 design coordinates (TP-2969)
morphgeom/                  the code (Python package)
  naca.py                   a = 1.0 mean line, 6-series assembly, NACA 4-digit (for tests), thickness tables
  surfaces.py               data types, camber + thickness -> surfaces, point loops, .dat reader and splitter
  morph.py                  v1 camber deflection, v2 constant_length, apply (picks the method), equivalent angles
  checks.py                 wing box, skin length, strain, joins, sanity
  aero.py                   rough lift change (thin-airfoil estimate)
  deflection.py             deflection range search, wing box size, one-setting-at-a-time study
  export.py                 .dat, CAD XYZ (mm), summary row / CSV, printed report
  plots.py                  overlay, camber lines, curvature, baseline check, deflection range, spar study, variable study
  pipeline.py               load_baseline, make_shape, write_outputs
  config.py                 reads config/default.yaml
  step.py                   STEP files (2D face, 3D solid) and a read-back check; needs pythonocc-core
scripts/                    make_baseline.py, make_shape.py, sweep.py, make_step.py, deflection_range.py, spar_study.py,
                            variable_study.py
tests/                      27 pytest tests (the 7 from the handoff, 5 helpers, 5 for v2, 3 for SC(2)-0714 / case files,
                            1 for STEP, 4 for the hinge, lift estimate, deflection range and spar study, 2 for the variable study)
outputs/                    generated files, one folder per airfoil:
  sc20714/                  SC(2)-0714 (default): baseline, shapes, sweep, range, spar_study, variable_study
  sc20714_tunnel/           SC(2)-0714 at 152.4 mm, including the STEP files
  naca641412/               NACA 64(1)-412: baseline, shapes, sweep, range, spar_study, variable_study
```

## Sources

- Thickness form (default): Abbott & von Doenhoff, *Theory of Wing Sections*, App. I (same as NACA Report 824), via the UIUC Airfoil Database file [n64012.dat](https://m-selig.ae.illinois.edu/ads/coord/n64012.dat).
- Reference 64(1)-412, 64(1)-012, 4412 and a = 1.0 tables: R. L. Carmichael, *Recomputed Tables of Coordinates of NACA Airfoil Profiles, Meanlines, and Sections*, Public Domain Aeronautical Software, [avd.pdf](https://www.pdas.com/avd.pdf) (pages 45, 104, 116, 147). Neither UIUC nor Airfoil Tools has a 64(1)-412 file.
- Ladson et al., NASA TM-4741, [program for NACA airfoil ordinates](https://www.pdas.com/refs/tm4741.pdf) (the method behind the PDAS tables).
- SC(2)-0714 coordinates: NASA TP-2969 (Harris, *NASA Supercritical Airfoils: A Matrix of Family-Related Airfoils*), via the UIUC Airfoil Database file [sc20714.dat](https://m-selig.ae.illinois.edu/ads/coord/sc20714.dat).
- SC(2)-0714 wind-tunnel data: Jenkins, Hill & Ray, NASA TM-4044, [Langley 0.3-m Transonic Cryogenic Tunnel](https://ntrs.nasa.gov/api/citations/19880015196/downloads/19880015196.pdf) (Mach 0.6–0.76); the same data corrected for sidewall boundary-layer effects: [NTRS 19950010493](https://ntrs.nasa.gov/api/citations/19950010493/downloads/19950010493.pdf).
- Morphing background and constant-arc-length LE: see the handoff's source list.
