Shapes for CFD: sc20714 baseline + 20 morphed trailing edges
Made by scripts/cfd_shapes.py from the cfd_set settings (config/default.yaml). Baseline: coordinates from sc20714.dat.

FILES
  <name>.dat       one shape each. Line 1 is the name, then one x y pair per line.
                   Chord = 1, 6 decimals, 200 points per surface (cosine spacing, bunched at the LE and TE).
                   Order: TE -> upper surface -> LE -> lower surface -> TE (Selig).
  cfd_cases.csv    every case: group, hinge, shape exponent n, flap angle, tip drop, rough lift estimate
  cfd_shapes.png   every shape drawn against the baseline

HOW THE SHAPES ARE SET UP
  - Only the trailing edge changes, behind the hinge ("xh" in the file name, x/c). Everything ahead
    of the hinge is the same as the baseline. The rear spar is at the hinge.
  - The flap angle ("flap+5deg") is the equivalent flap angle: a plain flap hinged at the same point
    with the same trailing-edge drop would sit at this angle. + = trailing edge down (more lift), - = up.
  - n is the shape exponent: how the bend is spread between the hinge and the trailing edge (the drop
    grows as distance^n). Higher n keeps the skin straighter near the hinge and bends it more near the tip.
  - The skins bend but keep their length (the lower skin slides a few mm into the trailing edge).
  - The shapes are NOT re-normalised. The baseline chord line stays on the x-axis and a drooped
    trailing edge simply sits lower. Run every shape at the same angle of attack, measured from the
    x-axis, so they are all compared at the same wing attitude.
  - Blunt trailing edge: 0.70% chord thick (7.0 mm at the placeholder 1 m chord).
    The first and last points are the top and bottom of the base; they are not joined.

FLIGHT CONDITION (config/default.yaml, flight:)
  Altitude 38,000 ft. Velocity 236.08 m/s, speed of sound 295.1 m/s: Mach 0.800.
  Temperature 216.7 K (from the speed of sound). Pressure 20,646 Pa. Density 0.3021 kg/m3 (as given).
  Reynolds number 5.02 million per metre of chord (5.02 million at the placeholder 1 m chord).
  Dynamic pressure 8,419 Pa.

STILL OPEN (question numbers from the README)
  - Q18 wing sweep: the section sees about Mach 0.80 x cos(sweep). Run at that Mach if the wing is swept.
  - Q17 density: the pressure and speed of sound give 0.332 kg/m3, not 0.3021. That would change the Reynolds number by +10%.
  - Q4 chord: 1 m is a placeholder, and the Reynolds number scales with it.
  - Q5 point format and how to close or mesh the trailing edge.
  - Angle of attack: to agree with aero (the cruise angle of attack, or a few angles around it).

PLEASE SEND BACK, for every shape and angle of attack
  - lift, drag and pitching moment coefficients (cl, cd, cm about the quarter chord)
  - the pressure distribution, where the shock sits, and whether the flow separates

THE GROUPS
  angle sweep (cases 2-9): hinge 0.650, n 3, flap angles -15, -10, -5, -2, +2, +5, +10, +13 deg
    to show the real lift per degree at cruise, and how far the trailing edge can go before shocks or separation stop it helping
  hinge (cases 10-15): hinge 0.550, 0.750, 0.850, n 3, flap angles -5, +5 deg
    to show how much the hinge position changes the lift (compare with hinge 0.65 in the angle sweep)
  shape exponent (cases 16-21): hinge 0.650, n 2.8, 3.5, 4, flap angles -5, +5 deg
    to show how much the shape exponent n changes the lift (compare with n 3 in the angle sweep)

The lift numbers in cfd_cases.csv are a rough low-speed estimate (thin-airfoil theory). They are not
reliable at this Mach number: replacing them is the point of these runs.
