"""Tests for the morphing shape generator. Run from the project folder:  pytest

Tests 1-7 are the ones listed in the handoff; the rest cover helpers and v2.
"""
import copy
from types import SimpleNamespace

import numpy as np
import pytest

from morphgeom import make_shape, write_outputs
from morphgeom.checks import bending_strain, curvature, joint_jumps, run_all, self_intersections
from morphgeom.export import shape_name, write_selig
from morphgeom.morph import apply, camber_delta, deflect_camber
from morphgeom.naca import camber_4digit, camber_a10, naca6_section, thickness_4digit
from morphgeom.surfaces import build, loop, read_dat, split_dat, stations


def morph(delta_LE, delta_TE, n, method="camber"):
    return SimpleNamespace(x_LE=0.15, x_h=0.65, delta_LE=delta_LE, delta_TE=delta_TE, n=n,
                           method=method, n_modes=4, le_smoothing=0.01, le_full_droop=True)


def with_morph(cfg, mp):
    c = copy.deepcopy(cfg)
    c.morph = mp
    return c


def surf(S):
    return np.column_stack([S.xu, S.yu, S.xl, S.yl])


# 1. Surface builder: NACA 0012 thickness + 4412 mean line = published NACA 4412
def test_builder_matches_published_naca4412(root, xref):
    S = build(camber_4digit(0.04, 0.4), thickness_4digit(0.12), xref)
    ref = read_dat(root / "data" / "naca4412_pdas_reference.dat")
    # PDAS gives 4 decimals in percent chord: rounding is at most 5e-7 per coordinate
    np.testing.assert_allclose(loop(S), ref.xy, rtol=0, atol=1e-6)


# 2. Zero morph returns the baseline exactly
def test_zero_morph_is_baseline(af):
    x = stations(200, "cosine")
    B = build(af.cam, af.thk, x)
    for n in (2, 3):
        M = build(deflect_camber(af.cam, morph(0, 0, n)), af.thk, x)
        assert np.array_equal(surf(M), surf(B))            # bit for bit


# 3. Wing box unchanged to machine precision
def test_wing_box_unchanged(af, cfg):
    mp = morph(0.01, 0.03, 3)
    x = stations(2001, "cosine", [mp.x_LE, mp.x_h])
    B = build(af.cam, af.thk, x)
    M = build(deflect_camber(af.cam, mp), af.thk, x)
    box = (x >= mp.x_LE) & (x <= mp.x_h)
    assert np.array_equal(surf(M)[box], surf(B)[box])
    assert np.max(np.abs(M.yu[~box] - B.yu[~box])) > 1e-3      # and something outside the box did move
    R = run_all(af, apply(af, mp, cfg), cfg)
    assert R["wingbox_station_max"] == 0


# 4. Joints: dy = 0 at the spars; slope jump ~0 for n >= 2; curvature jump ~0 for n >= 3
def test_joints(af):
    for n in (1, 2, 3, 4):
        mp = morph(0.01, 0.02, n)
        dy, _ = camber_delta(np.array([mp.x_LE, mp.x_h]), mp)
        assert np.array_equal(dy, [0.0, 0.0])
        cam = deflect_camber(af.cam, mp)
        shape = lambda x, cam=cam: build(cam, af.thk, x)  # noqa: E731
        Jh, Jle = joint_jumps(shape, mp.x_h), joint_jumps(shape, mp.x_LE)
        slope = max(Jh["slope_deg"][0], Jle["slope_deg"][0])      # camber line
        curv = max(Jh["curv"][0], Jle["curv"][0])
        if n == 1:
            assert slope > 1                                       # a kink: the check must see it
        else:
            assert slope < 1e-5
        if n == 2:
            # camber curvature jump at the rear spar = 2 delta / (1 - x_h)^2
            assert Jh["curv"][0] == pytest.approx(2 * mp.delta_TE / (1 - mp.x_h) ** 2, rel=1e-3)
        elif n >= 3:
            assert curv < 1e-5


# 5. Strain worked example: 2.000 mm skin, flat -> R = 0.500 m gives 0.002
def test_strain_worked_example():
    assert bending_strain(2.000, 1 / 0.500) == pytest.approx(0.002, abs=1e-15)


# 6. Baseline 64(1)-412 matches published coordinates
def test_baseline_matches_published(root, af, xref):
    ref = read_dat(root / "data" / "naca64_1-412_pdas_reference.dat")
    # default build (Abbott & von Doenhoff thickness): within 0.01% chord.
    # The two thickness tables differ by up to 0.006% chord, so this is the
    # real agreement, not a loose tolerance.
    S = build(af.cam, af.thk, xref)
    assert np.max(np.linalg.norm(loop(S) - ref.xy, axis=1)) < 1e-4
    # same thickness table PDAS used: only their 4-decimal rounding is left
    afp = naca6_section(root / "data" / "naca64_1-012_thickness_pdas.dat", 0.4, "p")
    np.testing.assert_allclose(loop(build(afp.cam, afp.thk, xref)), ref.xy, rtol=0, atol=1e-6)
    # mean line formula vs the PDAS a = 1.0 table (cl_i = 1.0)
    T = np.loadtxt(root / "data" / "meanline_a1.0_cl1.0_pdas.dat", comments="#")
    cam = camber_a10(1.0)
    np.testing.assert_allclose(cam.yc(T[:, 0]), T[:, 1], rtol=0, atol=1e-6)
    inner = (T[:, 0] > 0) & (T[:, 0] < 1)                 # slope is infinite at the ends (PDAS prints 0)
    np.testing.assert_allclose(cam.dyc(T[inner, 0]), T[inner, 2], rtol=0, atol=1e-4)


# 7. Output files re-parse in Selig order
def test_output_format(cfg, af, tmp_path):
    c = copy.deepcopy(cfg)
    c.output.plots = False
    c.structure.chord_m = 2.5
    res = make_shape(c, af)
    files = write_outputs(res, c, tmp_path, make_plots=False)

    d = read_dat(files["dat"])
    xy = loop(res.shape)
    N = c.points.n_per_surface
    assert d.format == "selig"
    assert d.name == res.name
    assert d.xy.shape == (2 * N - 1, 2)
    np.testing.assert_allclose(d.xy, xy, rtol=0, atol=5e-7)          # 6 decimals
    assert np.array_equal(d.xy[0], d.xy[-1])                         # closed TE: starts and ends at the TE
    assert abs(int(np.argmin(d.xy[:, 0])) - (N - 1)) < 20            # LE sits in the middle
    assert d.xy[:N, 1].mean() > d.xy[N - 1:, 1].mean()               # upper first, then lower
    np.testing.assert_allclose(d.xy[N - 1], [0, 0], atol=5e-7)       # station x = 0 point

    cad = np.loadtxt(files["cad"])
    np.testing.assert_allclose(cad, np.column_stack([xy * 2500, np.zeros(2 * N - 1)]), rtol=0, atol=5e-7)


# Extras
def test_file_name_matches_handoff_example():
    assert shape_name("naca641412", morph(0, 0.02, 3)) == "naca641412_TE+0.020_xh0.650_n3"
    assert shape_name("naca641412", morph(0, 0, 3)) == "naca641412_baseline"


def test_curvature_of_circle():
    t = np.linspace(0, np.pi, 2001)
    xy = 0.5 * np.column_stack([np.cos(t), np.sin(t)])         # counter-clockwise, R = 0.5
    np.testing.assert_allclose(curvature(xy)[1:-1], 2.0, rtol=1e-6)
    np.testing.assert_allclose(curvature(xy[::-1])[1:-1], -2.0, rtol=1e-6)   # clockwise: negative


def test_self_intersections(af):
    assert self_intersections(np.array([[0, 0], [1, 1], [1, 0], [0, 1], [0, 0]], float)) == 1   # bow tie
    assert self_intersections(np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]], float)) == 0   # square
    assert self_intersections(loop(build(af.cam, af.thk, stations(400, "cosine")))) == 0


def test_split_dat_rebuilds_file(root, tmp_path):
    a = split_dat(root / "data" / "naca64_1-412_pdas_reference.dat", "ref")
    assert a.split_error < 1e-5
    # symmetric airfoil: camber must come out flat
    S = build(camber_4digit(0, 0), thickness_4digit(0.12), stations(200, "cosine"))
    f = tmp_path / "naca0012.dat"
    write_selig(f, loop(S), "naca0012")
    a0 = split_dat(f, "n0012")
    assert np.max(np.abs(a0.cam.yc(np.linspace(0, 1, 101)))) < 1e-6
    assert a0.split_error < 1e-5


def test_checks_on_default_shape(cfg, af):
    # Sanity of the check values themselves on the default TE droop
    R = make_shape(cfg, af).checks
    assert R["arc_te_upper_pct"] > R["arc_te_lower_pct"]          # outer skin stretches more
    assert abs(R["arc_le_upper_pct"]) < 1e-12                     # LE untouched
    assert abs(R["strain_le_upper_pct"]) < 1e-9
    assert R["thickness_pass"] and R["intersect_pass"] and R["te_pass"] and R["le_pass"] and R["wingbox_pass"]
    # TE strain is close to the camber-line estimate (t/2) * 6 delta / (1 - x_h)^2 at the tip
    est = 100 * bending_strain(cfg.structure.t_skin_mm,
                               6 * cfg.morph.delta_TE / (1 - cfg.morph.x_h) ** 2 / cfg.structure.chord_m)
    assert 0.5 * est < R["strain_te_upper_pct"] < est


# v2: constant-length morph
def test_constlen_zero_morph_is_baseline(af, cfg):
    x = stations(200, "cosine")
    B = build(af.cam, af.thk, x)
    M = apply(af, morph(0, 0, 3, "constant_length"), cfg).eval(x)
    assert np.array_equal(surf(M), surf(B))                       # bit for bit


def test_constlen_keeps_skin_lengths(af, cfg):
    mp = morph(0.01, 0.02, 3, "constant_length")
    R = run_all(af, apply(af, mp, cfg), cfg)
    # every skin piece keeps its length, except the lower TE skin, which slides in
    # (1e-5 %: the checks sample the skin more coarsely than v2 builds it)
    for f in ("arc_te_upper_pct", "arc_le_upper_pct", "arc_le_lower_pct", "arc_le_wrap_pct"):
        assert abs(R[f]) < 1e-5, f
    assert R["te_lower_slide_mm"] > 1                            # a 2% droop pulls the lower skin in
    assert R["wingbox_station_max"] == 0                          # box untouched
    assert R["self_intersections"] == 0
    assert R["te_gap_change_mm"] == pytest.approx(0, abs=1e-9)   # TE still closed
    D = R["dense"]
    assert D["morph"].yu[-1] == pytest.approx(-mp.delta_TE, abs=1e-12)   # TE dropped by delta_TE
    assert D["morph"].xu[-1] < 1                                          # ... and moved forward
    # ribs keep their length behind the rear spar: thickness unchanged there
    te = D["x"] > mp.x_h
    tB = np.hypot(D["base"].xu - D["base"].xl, D["base"].yu - D["base"].yl)
    tM = np.hypot(D["morph"].xu - D["morph"].xl, D["morph"].yu - D["morph"].yl)
    assert np.max(np.abs(tM[te] - tB[te])) < 1e-12


def test_constlen_joints(af, cfg):
    # The n rule holds on the skins v2 is built from: slope matches for
    # n >= 2, curvature matches for n >= 3. v2 is a chain of short segments,
    # so "zero" here means the size of its discretisation error.
    for n in (2, 3):
        mr = apply(af, morph(0.01, 0.02, n, "constant_length"), cfg)
        Jh, Jle = joint_jumps(mr.eval, 0.65), joint_jumps(mr.eval, 0.15)
        assert max(Jh["slope_deg"][1], *Jle["slope_deg"][1:]) < 1e-3
        if n == 2:
            assert Jh["curv"][1] > 0.1                            # the check sees the curvature step
        else:
            assert max(Jh["curv"][1], *Jle["curv"][1:]) < cfg.checks.tol_curv_per_m


def test_constlen_rejects_kink(af, cfg):
    # n = 1 would put a kink in the skin at the spar
    with pytest.raises(ValueError):
        apply(af, morph(0, 0.02, 1, "constant_length"), cfg)


def test_constlen_reaches_full_droop(af, cfg):
    # The v2 nose drops exactly as far as the v1 nose (at 1% chord), keeps its
    # skin length, and bends with less strain than v1.
    for dl in (-0.01, 0.01):
        R1 = make_shape(with_morph(cfg, morph(dl, 0, 3)), af).checks
        R2 = make_shape(with_morph(cfg, morph(dl, 0, 3, "constant_length")), af).checks
        assert R2["nose_drop_mm"] == pytest.approx(R1["nose_drop_mm"], abs=1e-3)
        assert abs(R2["arc_le_wrap_pct"]) < 1e-5
        assert max(R2["strain_le_upper_pct"], R2["strain_le_lower_pct"]) < \
            max(R1["strain_le_upper_pct"], R1["strain_le_lower_pct"])
        assert R2["self_intersections"] == 0


# SC(2)-0714 (loaded .dat, blunt TE) and case files
@pytest.fixture(scope="module")
def sc_cfg():
    from morphgeom import load_config
    return load_config("config/sc20714.yaml")


def test_case_file_only_changes_what_it_lists(sc_cfg, cfg, tmp_path):
    from morphgeom import load_config
    assert sc_cfg.baseline.source == "dat" and sc_cfg.baseline.name == "sc20714"
    assert load_config().baseline.name == "sc20714"         # SC(2)-0714 is the default
    assert cfg.baseline.source == "naca6" and cfg.output.dir.endswith("outputs/naca641412")
    assert sc_cfg.morph.x_h == cfg.morph.x_h and sc_cfg.checks.n_points == cfg.checks.n_points
    bad = tmp_path / "bad.yaml"
    bad.write_text("morph:\n  x_hinge: 0.7\n")          # misspelt setting
    with pytest.raises(KeyError):
        load_config(bad)


def test_sc20714_loads_smooth_and_close_to_file(sc_cfg):
    from morphgeom import load_baseline
    af = load_baseline(sc_cfg)
    assert af.frame == "file frame kept"                 # NASA frame: angles of attack stay comparable
    assert af.split_error < 1e-4                         # within the file's rounding
    assert af.split_rms < 3e-5
    assert 2 * af.thk.te_half_thickness == pytest.approx(0.007, abs=1e-9)   # blunt TE kept exactly
    assert 0.025 < af.thk.le_radius_implied < 0.035
    # smoothed: no curvature sign changes on the nose skin of the baseline
    R = make_shape(with_morph(sc_cfg, morph(0, 0.0001, 3)), af).checks
    assert R["le_inflections_base"] == 0


def test_sc20714_v2_trailing_edge(sc_cfg):
    R = make_shape(with_morph(sc_cfg, morph(0, 0.02, 3, "constant_length"))).checks
    assert R["all_pass"]
    assert abs(R["arc_te_upper_pct"]) < 1e-5
    assert R["te_gap_change_mm"] < 1e-6                  # blunt TE gap unchanged (it turns with the TE)
    assert R["self_intersections"] == 0


# STEP files (skipped if pythonocc-core is not installed)
def test_step_face_and_solid_match_the_points(af, sc_cfg, tmp_path):
    step = pytest.importorskip("morphgeom.step")
    from morphgeom import load_baseline
    chord, span = 152.4, 100.0
    for a in (af, load_baseline(sc_cfg)):                # closed TE and blunt TE
        xy = loop(build(a.cam, a.thk, stations(200, "cosine"))) * chord
        area = 0.5 * abs(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], -1)))
        step.write_step(tmp_path / "face.step", xy)
        step.write_step(tmp_path / "solid.step", xy, span_mm=span)
        face, solid = step.check_step(tmp_path / "face.step"), step.check_step(tmp_path / "solid.step")
        assert face["valid"] and face["solids"] == 0 and face["faces"] == 1
        assert face["face_areas_mm2"][0] == pytest.approx(area, rel=1e-3)
        assert solid["valid"] and solid["solids"] == 1
        x0, y0, z0, x1, y1, z1 = solid["bbox_mm"]
        # the box is the points' own extent (the 64(1)-412 nose sits 0.035 mm ahead of x = 0)
        assert (x0, x1) == (pytest.approx(xy[:, 0].min(), abs=1e-3), pytest.approx(xy[:, 0].max(), abs=1e-3))
        assert (y0, y1) == (pytest.approx(xy[:, 1].min(), abs=1e-3), pytest.approx(xy[:, 1].max(), abs=1e-3))
        assert (z0, z1) == (pytest.approx(-span / 2, abs=1e-5), pytest.approx(span / 2, abs=1e-5))
        assert solid["volume_mm3"] == pytest.approx(area * span, rel=1e-3)


# Hinge behind the rear spar, deflection range, spar study, lift estimate
def test_hinge_behind_rear_spar(af, cfg):
    for method in ("camber", "constant_length"):
        mp = morph(0, 0.02, 3, method)
        mp.x_rear_spar, mp.x_h = 0.55, 0.70
        res = make_shape(with_morph(cfg, mp), af)
        R = res.checks
        assert R["wingbox_pass"] and R["slope_join_pass"] and R["curv_join_pass"]
        D = R["dense"]
        fixed = (D["x"] >= 0.15) & (D["x"] <= 0.70)          # wing box AND the skin up to the hinge
        assert np.array_equal(surf(D["morph"])[fixed], surf(D["base"])[fixed])
        assert "_xh0.700_rs0.550_" in res.name
    assert R["all_pass"]                                      # v2
    mp.x_rear_spar = 0.75                                     # spar behind the hinge: not allowed
    with pytest.raises(ValueError, match="hinge"):
        make_shape(with_morph(cfg, mp), af)


def test_lift_estimate(af):
    import math
    from morphgeom.aero import lift_change, zero_lift_cl
    # the a = 1.0 mean line with design cl 0.4 gives cl = 0.4 at zero angle of attack
    assert zero_lift_cl(build(af.cam, af.thk, stations(2001, "cosine"))) == pytest.approx(0.4, abs=1e-3)
    # n = 1 is a straight flap: textbook thin-airfoil flap lift
    mp = morph(0, 0.01, 1)
    mp.x_h = 0.7
    x = stations(2001, "cosine", [mp.x_LE, mp.x_h])
    dcl = lift_change(build(af.cam, af.thk, x), build(deflect_camber(af.cam, mp), af.thk, x))
    t = math.acos(1 - 2 * mp.x_h)
    assert dcl == pytest.approx(2 * (0.01 / 0.3) * (math.pi - t + math.sin(t)), rel=1e-9)


def test_deflection_range(af, cfg):
    from morphgeom import copy_config
    from morphgeom.deflection import find_range, try_angle
    c = copy_config(cfg)
    c.morph.method = "constant_length"
    rng = find_range(c, af, "TE")
    for lim in rng.values():
        assert lim.stopped_by == ["strain"]
        strain = lim.res.checks["strain_max_pct"]
        assert 0.98 * c.checks.tol_strain_pct < strain <= c.checks.tol_strain_pct
        assert try_angle(c, af, "TE", lim.angle_deg + np.sign(lim.angle_deg) * c.range.tol_deg)[1]   # just past: fails
    assert rng["droop"].angle_deg > 5 and rng["reflex"].angle_deg < -5
    c.morph.method = "camber"                                 # v1 stretches its skins
    assert find_range(c, af, "TE")["droop"].stopped_by == ["skin_length"]
    c.range.ignore = ["skin_lenght"]                          # misspelt check name
    with pytest.raises(KeyError):
        find_range(c, af, "TE")


def test_spar_further_back_reaches_less(af, cfg):
    from morphgeom import copy_config
    from morphgeom.deflection import find_range, wing_box
    out = []
    for xs in (0.60, 0.75):
        c = copy_config(cfg)
        c.morph.method, c.morph.x_h = "constant_length", xs
        out.append((find_range(c, af, "TE")["droop"].angle_deg, wing_box(af, c.morph.x_LE, xs)))
    (d1, b1), (d2, b2) = out
    assert d2 < d1                                            # less skin to bend: less droop
    assert b2["rear_height"] < b1["rear_height"] and b2["area_share"] > b1["area_share"]


def test_strain_settings_act_as_one_number(af, cfg):
    """Chord, skin thickness and strain limit reach the range only through
    skin / (2 x strain limit) / chord: doubling any one the right way gives
    the same limit."""
    from morphgeom import copy_config
    from morphgeom.deflection import bend_radius_limit, find_limit
    c = copy_config(cfg)
    c.morph.method = "constant_length"
    c.range.step_deg, c.range.tol_deg = 2.0, 0.05
    out = []
    for path, v in (("structure.chord_m", 2 * c.structure.chord_m),
                    ("checks.tol_strain_pct", 2 * c.checks.tol_strain_pct),
                    ("structure.t_skin_mm", c.structure.t_skin_mm / 2)):
        ci = copy_config(c)
        sec, key = path.split(".")
        setattr(getattr(ci, sec), key, v)
        assert bend_radius_limit(ci) == pytest.approx(bend_radius_limit(c) / 2)
        out.append(find_limit(ci, af, "TE", +1).angle_deg)
    assert max(out) - min(out) < 1e-9                    # the very same shapes pass and fail
    assert out[0] > find_limit(c, af, "TE", +1).angle_deg + 1   # and it really is further


def test_vary_one(af, cfg):
    from morphgeom import copy_config
    from morphgeom.deflection import get_setting, set_setting, vary_one
    c = copy_config(cfg)
    c.morph.method = "constant_length"
    c.range.step_deg, c.range.tol_deg = 2.0, 0.1
    cache = {}
    rows = vary_one(c, af, "TE", {"setting": "morph.n", "values": [3, 4.0], "label": "n"}, cache)
    assert [r["value"] for r in rows] == ["3", "4"] and [r["is_default"] for r in rows] == ["yes", "no"]
    assert rows[1]["droop_max_deg"] < rows[0]["droop_max_deg"]   # n = 4 bends harder near the tip
    assert len(cache) == 2
    again = vary_one(c, af, "TE", {"setting": "morph.n", "values": [3.0]}, cache)
    assert len(cache) == 2 and again[0]["droop_max_deg"] == rows[0]["droop_max_deg"]   # 3.0 = 3: reused
    # the rear spar ahead of a fixed hinge does not change the TE range
    held = vary_one(c, af, "TE", {"setting": "morph.x_rear_spar", "values": [0.55]}, cache)
    assert held[0]["droop_max_deg"] == pytest.approx(rows[0]["droop_max_deg"], abs=1e-9)
    set_setting(c, "morph.n", 5)
    assert get_setting(c, "morph.n") == 5
    with pytest.raises(KeyError):
        vary_one(c, af, "TE", {"setting": "morph.x_hinge", "values": [0.7]})


def test_flight_condition():
    from morphgeom.aero import flight_condition, lift_per_span
    fl = SimpleNamespace(speed_of_sound_m_s=295.1, velocity_m_s=236.08, density_kg_m3=0.3021, pressure_pa=20646)
    fc = flight_condition(fl)
    assert fc["mach"] == pytest.approx(0.8, abs=1e-4)
    assert fc["q_pa"] == pytest.approx(0.5 * 0.3021 * 236.08 ** 2)
    assert fc["prandtl_glauert"] == pytest.approx(1 / 0.6, rel=1e-3)
    # standard atmosphere at 38,000 ft: 216.65 K and 0.3320 kg/m3 (so 0.3021 does not match)
    assert fc["temperature_k"] == pytest.approx(216.65, abs=0.1)
    assert fc["density_from_p_t"] == pytest.approx(0.3320, abs=5e-4)
    assert lift_per_span(0.1, 8000.0, 1.5) == pytest.approx(1200.0)


def test_flap_point(af, cfg):
    from morphgeom import copy_config
    from morphgeom.deflection import find_range, flap_point
    c = copy_config(cfg)
    c.morph.method = "constant_length"
    c.range.step_deg, c.range.tol_deg = 2.0, 0.1
    rows, lim = flap_point(c, 0.7, 4, [-2.0, 2.0, 40.0], af)
    assert [r["angle_deg"] for r in rows] == [-2.0, 2.0, 40.0]
    assert rows[0]["passes"] and rows[1]["passes"] and not rows[2]["passes"] and rows[2]["failed"]
    assert rows[1]["dcl_estimate"] > 0 > rows[0]["dcl_estimate"]          # droop adds lift
    assert rows[1]["tip_mm"] == pytest.approx(0.3 * np.tan(np.radians(2)) * 1000 * c.structure.chord_m)
    c.morph.x_h, c.morph.n = 0.7, 4                                       # same as find_range directly
    rng = find_range(c, af, "TE")
    assert lim["droop_max_deg"] == rng["droop"].angle_deg and lim["reflex_max_deg"] == rng["reflex"].angle_deg
    assert lim["range_deg"] == pytest.approx(lim["droop_max_deg"] - lim["reflex_max_deg"])
