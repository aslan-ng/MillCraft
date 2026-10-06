from collections import Counter
import math

import cadquery as cq
import numpy as np
from PIL import Image
import pytest

from millcraft.cli import generate, parser
from millcraft.machining import envelope_depth, export_obj, machine_plate
from millcraft.sampling import MM_PER_INCH, Parameters, load_grayscale, sample_image
from millcraft.splines import interpolate_grid, interpolate_row


def small_plate():
    return Parameters(plate_width=0.5, plate_height=0.5, plate_thickness=0.1,
                      artwork_width=0.3, artwork_height=0.3, max_depth=0.02)


def flat_row(y=0, depth=-0.508):
    return interpolate_row([-3, 3], y, [depth, depth], 0.508)


def test_constant_depth_capsule_volume_and_round_ends():
    p = small_plate()
    plate, report = machine_plate([flat_row()], p)
    r = p.tool_diameter * MM_PER_INCH / 2
    polygon_area = 32 * r**2 * math.sin(2 * math.pi / 32) / 2
    expected = (6 * 2 * r + polygon_area) * 0.508
    assert report['removed_volume_mm3'] == pytest.approx(expected, abs=1e-5)
    assert plate.isValid() and len(plate.Solids()) == 1
    # Flat floors, rounded end footprints, and untouched material between rows.
    assert not plate.isInside((0, 0, -0.2))
    assert not plate.isInside((3.5, 0, -0.2))
    assert plate.isInside((3.7, 0.7, -0.2))
    assert plate.isInside((0, r + 0.1, -0.2))
    assert plate.isInside((0, 0, -0.6))


def test_radius_envelope_uses_deeper_neighbor_and_analytic_valley():
    row = interpolate_row([-3, 3], 0, [-0.508, 0], 0.508)
    assert envelope_depth(row, np.array([0.0]), 0.5)[0] == pytest.approx(row.polynomial(-0.5))
    valley = interpolate_row([-3, 0, 3], 0, [0, -0.508, 0], 0.508)
    assert envelope_depth(valley, np.array([0.2]), 0.5)[0] == pytest.approx(-0.508)


def test_overlapping_grooves_remove_union_and_white_row_is_skipped():
    p = small_plate()
    single, one = machine_plate([flat_row()], p)
    overlap, two = machine_plate([flat_row(-0.4), flat_row(0.4)], p)
    assert overlap.isValid() and len(overlap.Solids()) == 1
    assert one['removed_volume_mm3'] < two['removed_volume_mm3'] < 2 * one['removed_volume_mm3']
    white, empty = machine_plate([flat_row(depth=0)], p)
    assert empty['active_grooves'] == 0
    assert empty['removed_volume_mm3'] == pytest.approx(0)
    assert white.Volume() > single.Volume()


def read_obj(path):
    vertices, faces = [], []
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields and fields[0] == 'v':
            vertices.append([float(value) for value in fields[1:]])
        elif fields and fields[0] == 'f':
            faces.append([int(value) - 1 for value in fields[1:]])
    return np.asarray(vertices), np.asarray(faces)


def test_obj_is_closed_oriented_and_matches_step_solid(tmp_path):
    plate, report = machine_plate([flat_row()], small_plate())
    step = tmp_path / 'plate.step'
    obj = tmp_path / 'plate.obj'
    cq.exporters.export(plate, str(step), exportType='STEP')
    imported = cq.importers.importStep(str(step)).solids().vals()
    assert len(imported) == 1 and imported[0].isValid()
    assert imported[0].Volume() == pytest.approx(report['final_volume_mm3'], abs=1e-5)
    export_obj(plate, obj)
    vertices, faces = read_obj(obj)
    counts = Counter(tuple(sorted((a, b))) for face in faces
                     for a, b in zip(face, np.roll(face, -1)))
    assert set(counts.values()) == {2}
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    mesh_volume = np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6
    assert mesh_volume == pytest.approx(plate.Volume(), abs=1e-4)
    bounds = plate.BoundingBox()
    np.testing.assert_allclose(vertices.min(axis=0), [bounds.xmin, bounds.ymin, bounds.zmin], atol=1e-6)
    np.testing.assert_allclose(vertices.max(axis=0), [bounds.xmax, bounds.ymax, bounds.zmax], atol=1e-6)


def test_default_cli_produces_one_step_and_matching_obj(tmp_path):
    source = tmp_path / 'black.png'
    Image.new('L', (8, 8), 0).save(source)
    output = tmp_path / 'result'
    args = parser().parse_args([str(source), '--out-dir', str(output),
                               '--plate-width', '0.5', '--plate-height', '0.5',
                               '--artwork-width', '0.3', '--artwork-height', '0.3',
                               '--line-spacing', '0.15', '--sample-spacing', '0.15'])
    report = generate(args)
    assert report['geometry'] == 'machined_plate_solid'
    assert report['machining']['removed_volume_mm3'] > 0
    assert len(list(output.glob('*.step'))) == 1
    assert (output / 'machined_plate.obj').exists()
    old_step = (output / 'machined_plate.step').read_bytes()
    args.cutter_segments = 7
    with pytest.raises(ValueError, match='Cutter segments'):
        generate(args)
    assert (output / 'machined_plate.step').read_bytes() == old_step


@pytest.mark.parametrize('min_depth', [0.0, 0.01, 0.035])
def test_white_image_uses_minimum_depth_in_finished_exports(tmp_path, min_depth):
    source = tmp_path / 'white.png'
    Image.new('L', (8, 8), 255).save(source)
    output = tmp_path / 'result'
    args = parser().parse_args([str(source), '--out-dir', str(output),
                               '--plate-width', '0.5', '--plate-height', '0.5',
                               '--artwork-width', '0.3', '--artwork-height', '0.3',
                               '--line-spacing', '0.15', '--sample-spacing', '0.15',
                               '--min-depth', str(min_depth)])
    report = generate(args)
    assert report['parameters_inch']['min_depth'] == min_depth
    assert report['validation']['z_max_mm'] == pytest.approx(-min_depth * MM_PER_INCH)
    assert report['validation']['minimum_depth_undershoot_mm'] <= 1e-7
    assert report['obj_mesh']['watertight']
    plate = cq.importers.importStep(str(output / 'machined_plate.step')).solids().val()
    assert plate.isValid()
    assert plate.isInside((0, 0, -min_depth * MM_PER_INCH - 0.02))
    if min_depth > 0:
        assert not plate.isInside((0, 0, -min_depth * MM_PER_INCH / 2))
        assert report['machining']['removed_volume_mm3'] > 0
    else:
        assert plate.isInside((0, 0, -0.02))
        assert report['machining']['removed_volume_mm3'] == pytest.approx(0)
    preview = (output / 'preview.svg').read_text()
    assert f'white = {min_depth:g} in depth' in preview
    assert 'nan' not in preview.lower()


def test_grazing_cut_on_supplied_image_is_a_valid_solid():
    from pathlib import Path

    source = Path(__file__).resolve().parents[1] / 'image.png'
    if not source.exists():
        pytest.skip('Source image is unavailable')
    p = Parameters(max_depth=0.20)
    rows = interpolate_grid(sample_image(load_grayscale(source), p), p.max_depth * MM_PER_INCH)
    plate, report = machine_plate([rows[2]], p)
    assert plate.isValid() and len(plate.Solids()) == 1
    assert report['removed_volume_mm3'] > 0


def test_switching_modes_retires_only_previous_generated_geometry(tmp_path):
    source = tmp_path / 'white.png'
    Image.new('L', (8, 8), 255).save(source)
    output = tmp_path / 'result'
    args = parser().parse_args([str(source), '--out-dir', str(output),
                               '--geometry', 'centerlines', '--with-plate',
                               '--line-spacing', '1.5', '--sample-spacing', '1.0'])
    generate(args)
    (output / 'user.step').write_text('user-owned file')
    args.geometry = 'machined'
    generate(args)
    assert not (output / 'centerlines.step').exists()
    assert not (output / 'plate_reference.step').exists()
    assert (output / 'user.step').read_text() == 'user-owned file'
    assert (output / 'machined_plate.step').exists()


@pytest.mark.parametrize('settings', [
    {'cut_sample_spacing': 0}, {'cut_sample_spacing': float('nan')},
    {'cutter_segments': 6}, {'cutter_segments': 130}, {'cut_sample_spacing': 1e-9},
])
def test_invalid_cut_resolution(settings):
    with pytest.raises(ValueError):
        machine_plate([flat_row()], small_plate(), **settings)
