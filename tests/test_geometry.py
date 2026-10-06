from dataclasses import replace
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GeomAbs import GeomAbs_BSplineCurve
from PIL import Image
import pytest

from millcraft.cli import generate, main, parser
from millcraft.sampling import MM_PER_INCH, Parameters, inclusive_axis, load_grayscale, sample_image
from millcraft.splines import SplineOvershootError, interpolate_grid, interpolate_row
from millcraft.step import export_centerlines, export_plate


def test_defaults_grid_mapping_and_orientation():
    p = Parameters()
    grid = sample_image(np.array([[0., 0.], [1., 1.]]), p)
    assert grid.z_mm.shape == (41, 121)
    np.testing.assert_allclose(np.diff(grid.x_mm), 0.025 * MM_PER_INCH)
    np.testing.assert_allclose(np.diff(grid.y_mm), 0.075 * MM_PER_INCH)
    np.testing.assert_allclose(grid.x_mm[[0, -1]], [-38.1, 38.1])
    np.testing.assert_allclose(grid.y_mm[[0, -1]], [-38.1, 38.1])
    np.testing.assert_allclose(grid.z_mm[0], 0)  # white bottom
    np.testing.assert_allclose(grid.z_mm[-1], -0.035 * MM_PER_INCH)  # black top
    np.testing.assert_allclose(grid.z_mm[20], -0.035 * MM_PER_INCH / 2)


def test_minimum_depth_maps_white_gray_and_black():
    p = Parameters(min_depth=0.01, max_depth=0.035)
    grid = sample_image(np.array([[1., 0.5, 0.], [1., 0.5, 0.]]), p, "stretch")
    np.testing.assert_allclose(grid.z_mm[:, 0], -0.01 * MM_PER_INCH)
    np.testing.assert_allclose(grid.z_mm[:, 60], -0.0225 * MM_PER_INCH)
    np.testing.assert_allclose(grid.z_mm[:, -1], -0.035 * MM_PER_INCH)
    rows = interpolate_grid(grid, p.max_depth * MM_PER_INCH,
                            min_depth_mm=p.min_depth * MM_PER_INCH)
    assert all(row.z_max_mm <= -p.min_depth * MM_PER_INCH + 1e-7 for row in rows)
    assert all(row.z_min_mm >= -p.max_depth * MM_PER_INCH - 1e-7 for row in rows)


def test_minimum_depth_applies_to_white_containment_margins():
    p = Parameters(min_depth=0.01)
    grid = sample_image(np.zeros((4, 2)), p)
    np.testing.assert_allclose(grid.z_mm[:, [0, -1]], -0.01 * MM_PER_INCH)
    np.testing.assert_allclose(grid.z_mm[:, 60], -p.max_depth * MM_PER_INCH)


def test_equal_depth_limits_produce_constant_depth():
    p = Parameters(min_depth=0.035, max_depth=0.035)
    grid = sample_image(np.array([[0., 1.], [0.5, 0.]]), p)
    np.testing.assert_allclose(grid.z_mm, -0.035 * MM_PER_INCH)


def test_bilinear_sampling_and_aspect_ratio():
    grid = sample_image(np.array([[0., 1.], [1., 1.]]), Parameters())
    assert grid.z_mm[20, 60] == pytest.approx(-0.035 * MM_PER_INCH * 0.25)
    narrow = np.zeros((4, 2))
    contained = sample_image(narrow, Parameters())
    np.testing.assert_allclose(contained.z_mm[:, 0], 0)
    np.testing.assert_allclose(contained.z_mm[:, 60], -0.035 * MM_PER_INCH)
    stretched = sample_image(narrow, Parameters(), "stretch")
    np.testing.assert_allclose(stretched.z_mm, -0.035 * MM_PER_INCH)


def test_transparency_is_white_and_blur_is_optional(tmp_path):
    image = tmp_path / "transparent.png"
    Image.new("RGBA", (2, 2), (0, 0, 0, 0)).save(image)
    np.testing.assert_array_equal(load_grayscale(image), np.ones((2, 2)))
    with pytest.raises(ValueError):
        load_grayscale(image, -1)


def test_nondivisible_spacing_preserves_step_and_endpoints():
    np.testing.assert_allclose(inclusive_axis(10, 3), [-5, -2, 1, 4, 5])
    np.testing.assert_allclose(inclusive_axis(10, 20), [-5, 5])


@pytest.mark.parametrize("changes", [
    {"line_spacing": 0}, {"sample_spacing": float("nan")}, {"max_depth": 0.25},
    {"artwork_width": 4}, {"tool_diameter": -1}, {"sample_spacing": 1e-9},
    {"min_depth": -0.01}, {"min_depth": float("nan")},
    {"min_depth": float("inf")}, {"min_depth": 0.036},
])
def test_invalid_parameters(changes):
    with pytest.raises(ValueError):
        replace(Parameters(), **changes)


def test_shape_preserving_spline_and_analytic_overshoot_rejection():
    x = np.array([0., 1., 2., 3.])
    z = np.array([0., -1., -1., 0.])
    bounded = interpolate_row(x, 0, z, 1, "pchip")
    np.testing.assert_allclose(bounded.polynomial(x), z)
    assert bounded.z_min_mm == pytest.approx(-1)
    assert bounded.z_max_mm == pytest.approx(0)
    assert bounded.local_overshoot_mm < 1e-12
    # Natural cubic dips below -1 BETWEEN two equally deep samples.
    with pytest.raises(SplineOvershootError, match="exceeds"):
        interpolate_row(x, 0, z, 1, "cubic")
    # Detect adjacent-sample overshoot even if it stays inside global limits.
    with pytest.raises(SplineOvershootError, match="adjacent"):
        interpolate_row(x, 0, z * 0.5 - 0.1, 1, "cubic")


def test_safe_natural_cubic_is_supported():
    row = interpolate_row([0, 1, 2, 3], 1, [0, -0.1, -0.2, -0.3], 1, "cubic")
    np.testing.assert_allclose(row.polynomial([0.5, 1.5, 2.5]), [-0.05, -0.15, -0.25])


def test_spline_rejects_depth_shallower_than_minimum():
    with pytest.raises(SplineOvershootError, match="exceeds"):
        interpolate_row([0, 1], 0, [-0.5, -0.2], 1, min_depth_mm=0.3)


def test_step_round_trip_entire_a_grid_and_plate(tmp_path):
    source = Path(__file__).resolve().parents[1] / "image.png"
    if not source.exists():
        pytest.skip("The user's source image is not in this checkout")
    parameters = Parameters()
    grid = sample_image(load_grayscale(source), parameters)
    rows = interpolate_grid(grid, parameters.max_depth * MM_PER_INCH)
    step_path = tmp_path / "a.step"
    export_centerlines(rows, step_path)
    imported = cq.importers.importStep(str(step_path))
    assert not imported.solids().vals()
    edges = sorted(imported.edges().vals(), key=lambda edge: edge.Center().y)
    assert len(edges) == 41
    assert "SI_UNIT(.MILLI.,.METRE.)" in step_path.read_text()
    for edge, row in zip(edges, rows):
        assert edge.isValid()
        curve = BRepAdaptor_Curve(edge.wrapped)
        assert curve.GetType() == GeomAbs_BSplineCurve
        assert curve.FirstParameter() == pytest.approx(-38.1)
        assert curve.LastParameter() == pytest.approx(38.1)
        # Round trip every sample plus ALL analytic extrema (not just dense probing).
        for x in np.unique(np.concatenate((grid.x_mm, row.extrema_x_mm))):
            point = curve.Value(float(x))
            np.testing.assert_allclose([point.X(), point.Y(), point.Z()],
                                       [x, row.y_mm, row.polynomial(x)], atol=1e-7, rtol=0)
    plate_path = tmp_path / "plate.step"
    export_plate(parameters, plate_path)
    solid = cq.importers.importStep(str(plate_path)).solids().val()
    box = solid.BoundingBox()
    np.testing.assert_allclose([box.xlen, box.ylen, box.zlen, box.zmin, box.zmax],
                              [88.9, 88.9, 6.35, -6.35, 0], atol=1e-7)
    assert solid.Volume() == pytest.approx(88.9 * 88.9 * 6.35)


def test_cli_outputs_and_failure_does_not_overwrite_step(tmp_path):
    source = tmp_path / "black.png"
    Image.new("L", (8, 8), 0).save(source)
    out = tmp_path / "results"
    args = parser().parse_args([str(source), "--out-dir", str(out), "--geometry", "centerlines", "--with-plate",
                               "--line-spacing", "1.5", "--sample-spacing", "1.0"])
    report = generate(args)
    assert report["grid"]["rows"] == 3
    assert report["grid"]["samples_per_row"] == 4
    assert report["validation"]["z_min_mm"] == pytest.approx(-0.889)
    assert all((out / name).is_file() for name in report["outputs"])
    old_step = (out / "centerlines.step").read_bytes()
    assert main([str(source), "--out-dir", str(out), "--max-depth", "1"]) == 1
    assert (out / "centerlines.step").read_bytes() == old_step
