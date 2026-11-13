naca64_1-412_pdas_reference.dat
-------------------------------
NACA 64(1)-412 section ordinates used to check the baseline build.

Source: R. L. Carmichael, "Recomputed Tables of Coordinates of NACA Airfoil Profiles,
Meanlines, and Sections" (based on NACA Report 824 and Abbott & von Doenhoff),
Public Domain Aeronautical Software, edition of 2022-05-17, https://www.pdas.com/avd.pdf,
page 147, "NACA Section 64-412", right-hand table (thickness added perpendicular to the
mean line, as in the book). PDAS drops the subscript: "64-412" here is the 64(1)-412.
Values were given in percent chord with 4 decimals; this file is in fraction of chord.

Why PDAS and not UIUC / Airfoil Tools: neither database has a 64(1)-412 file
(checked 2026-10-22).

Format: Lednicer (upper surface LE->TE, then lower surface LE->TE). 26 points each.

naca4412_pdas_reference.dat
---------------------------
NACA 4412 section, same source (avd.pdf page 116, "NACA Section 4412", right-hand table).
Standard 4-digit thickness (open trailing edge, yt(1) = 0.00126). Used by the tests to
check the surface builder against published numbers.

downloads/sc20714.dat
---------------------
NASA SC(2)-0714 supercritical airfoil, design coordinates from NASA TP-2969
(Harris, "NASA Supercritical Airfoils: A Matrix of Family-Related Airfoils", 1990),
as stored in the UIUC Airfoil Coordinates Database:
https://m-selig.ae.illinois.edu/ads/coord/sc20714.dat (downloaded 2026-10-27).
Lednicer format, 103 points per surface, blunt trailing edge (gap 0.0070 c),
TE midpoint at y = -0.013. Printed to 5 decimals; the points scatter like
4-decimal values, so the loader smooths to that level (README, "NASA SC(2)-0714").
Used by config/sc20714.yaml. Kept exactly as downloaded.
