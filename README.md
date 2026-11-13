# Morphing wing geometry

Shape generator for a morphing wing (ERAU geometry team). It takes a baseline airfoil,
bends the leading and trailing edges, checks each shape against structural limits,
and writes coordinate files and STEP files for SolidWorks.

**Start here: [morphing-geometry/README.md](morphing-geometry/README.md)** (set-up, commands, results, open questions).

| Folder | What is in it |
|---|---|
| `morphing-geometry/` | The Python tools, case files, tests and generated outputs |
| `cad-files/` | SolidWorks parts made from the generated shapes |

> The chord, skin thickness, strain limit, spar positions and needed deflection range
> are still **placeholders** until the structures and aero teams supply real values,
> so the results in `outputs/` show the method, not final design numbers.
