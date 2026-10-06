"""Single-spline row pockets and an optional cylindrical cutter-envelope model."""

from collections import Counter
import math
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCP.BRepGProp import BRepGProp
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepTools import BRepTools
from OCP.Geom import Geom_BSplineCurve
from OCP.GProp import GProp_GProps
from OCP.TColgp import TColgp_Array1OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.gp import gp_Pnt

from .sampling import MM_PER_INCH
from .step import plate_shape, row_edge

SURFACE_CLEARANCE_MM = 0.0001


def linear_edge(points):
    """A single C0 B-spline edge through a polyline, without depth overshoot."""
    count = len(points)
    poles = TColgp_Array1OfPnt(1, count)
    knots = TColStd_Array1OfReal(1, count)
    multiplicities = TColStd_Array1OfInteger(1, count)
    for i, point in enumerate(points, 1):
        poles.SetValue(i, gp_Pnt(*map(float, point)))
        knots.SetValue(i, (i - 1) / (count - 1))
        multiplicities.SetValue(i, 2 if i in (1, count) else 1)
    curve = Geom_BSplineCurve(poles, knots, multiplicities, 1, False)
    builder = BRepBuilderAPI_MakeEdge(curve)
    if not builder.IsDone():
        raise RuntimeError("Could not construct cutter envelope edge")
    return cq.Edge(builder.Edge())


def envelope_depth(row, x, half_width):
    """Exact spline minima within the cutter footprint at section coordinates."""
    curve = row.polynomial
    a, b = curve.x[[0, -1]]
    lo = np.clip(np.asarray(x) - half_width, a, b)
    hi = np.clip(np.asarray(x) + half_width, a, b)
    z = np.minimum(curve(lo), curve(hi))
    # Includes every analytic extremum, so narrow valleys cannot be missed
    # merely because they lie between source-image sampling knots.
    for t in row.extrema_x_mm:
        z = np.where((lo <= t) & (t <= hi), np.minimum(z, curve(t)), z)
    return z


def cutter_volume(row, radius_mm, spacing_mm, segments):
    a, b = row.polynomial.x[[0, -1]]
    count = math.ceil((b - a + 2 * radius_mm) / spacing_mm) + 1
    fraction = np.linspace(0, 1, count)
    wires = []
    # These sections trace an inscribed regular polygon around the circular
    # cutter. Every floor section shares the same normalized knot vector.
    for angle in np.linspace(-math.pi / 2, math.pi / 2, segments // 2 + 1):
        half_width = max(0.0, radius_mm * math.cos(angle))
        y = row.y_mm + radius_mm * math.sin(angle)
        x = a - half_width + fraction * (b - a + 2 * half_width)
        z = envelope_depth(row, x, half_width)
        # Avoid grazing/coincident faces at white image regions. Only the
        # uppermost micron is affected; deeper groove floors stay unchanged.
        z += SURFACE_CLEARANCE_MM * np.clip(1 + z / 0.001, 0, 1)
        points = np.column_stack((x, np.full(count, y), z))
        first, last = tuple(points[0]), tuple(points[-1])
        # All cutter sections extend above the plate, including white regions.
        top_first, top_last = (first[0], y, 1.0), (last[0], y, 1.0)
        wires.append(cq.Wire.assembleEdges([
            linear_edge(points),
            cq.Edge.makeLine(last, top_last),
            cq.Edge.makeLine(top_last, top_first),
            cq.Edge.makeLine(top_first, first),
        ]))
    tool = cq.Solid.makeLoft(wires, ruled=True)
    if not tool.isValid():
        raise RuntimeError(f"Invalid cutter volume at Y={row.y_mm:.6f} mm")
    return tool


def spline_pocket_volume(row, width_mm):
    """Extrude one exact row spline into a rectangular pocket with straight ends."""
    half_width = width_mm / 2
    floor = row_edge(row).translate((0, -half_width, 0))
    first, last = floor.startPoint(), floor.endPoint()
    top_start = cq.Vector(first.x, first.y, 1)
    top_finish = cq.Vector(last.x, last.y, 1)
    profile = cq.Wire.assembleEdges([
        floor, cq.Edge.makeLine(last, top_finish),
        cq.Edge.makeLine(top_finish, top_start), cq.Edge.makeLine(top_start, first),
    ])
    pocket = cq.Solid.extrudeLinear(profile, [], (0, width_mm, 0))
    if not pocket.isValid() or len(pocket.Solids()) != 1:
        raise RuntimeError(f"Invalid spline pocket at Y={row.y_mm:.6f} mm")
    return pocket


def solid_volume(shape):
    """Use adaptive integration across spline spans for accurate volume reports."""
    properties = GProp_GProps()
    error = BRepGProp.VolumePropertiesGK_s(shape.wrapped, properties, 1e-9, True, True)
    if error < 0 or not math.isfinite(properties.Mass()):
        raise RuntimeError("Could not calculate solid volume")
    return properties.Mass()


def machine_plate(rows, parameters, cut_sample_spacing=0.01, cutter_segments=32,
                  progress=None, cut_model="spline-pocket"):
    if cut_model not in ("spline-pocket", "cutter-envelope"):
        raise ValueError("Cut model must be 'spline-pocket' or 'cutter-envelope'")
    radius = parameters.tool_diameter * MM_PER_INCH / 2
    model_report = {
        "model": "single_spline_rectangular_pocket",
        "approximate": False,
        "floor_interpolation": "one exact row spline extruded uniformly across Y",
        "pocket_width_mm": radius * 2,
        "end_shape": "straight; at the spline X endpoints",
        "surface_clearance_mm": 0.0,
    }
    if cut_model == "cutter-envelope":
        if not math.isfinite(cut_sample_spacing) or cut_sample_spacing <= 0:
            raise ValueError("Cut sample spacing must be finite and positive (inches)")
        if not isinstance(cutter_segments, int) or not 8 <= cutter_segments <= 128 or cutter_segments % 4:
            raise ValueError("Cutter segments must be a multiple of 4 between 8 and 128")
        spacing = cut_sample_spacing * MM_PER_INCH
        section_samples = math.ceil((parameters.artwork_width * MM_PER_INCH + 2 * radius) / spacing) + 1
        if section_samples * (cutter_segments // 2 + 1) * len(rows) > 2_000_000:
            raise ValueError("Requested cutter discretization exceeds the 2,000,000 point limit")
        model_report = {
            "model": "vertical_cylindrical_flat_end_mill_envelope",
            "approximate": True,
            "cut_sample_spacing_mm": spacing,
            "cutter_segments": cutter_segments,
            "cutter_radial_chord_error_mm": radius * (1 - math.cos(math.pi / cutter_segments)),
            "floor_interpolation": "piecewise linear in X; ruled between transverse sections",
            "surface_clearance_mm": SURFACE_CLEARANCE_MM,
        }
    plate = plate_shape(parameters)
    original_volume = solid_volume(plate)
    threshold = SURFACE_CLEARANCE_MM if cut_model == "cutter-envelope" else 1e-7
    active = [row for row in rows if row.z_min_mm < -threshold]
    for i, row in enumerate(active, 1):
        if progress:
            progress(f"Cutting groove {i}/{len(active)}...")
        pocket = (spline_pocket_volume(row, radius * 2) if cut_model == "spline-pocket"
                  else cutter_volume(row, radius, spacing, cutter_segments))
        plate = plate.cut(pocket)
    if not plate.isValid() or len(plate.Solids()) != 1:
        raise RuntimeError("Machining did not produce one valid plate solid")
    plate = plate.Solids()[0]
    bounds = plate.BoundingBox()
    thickness = parameters.plate_thickness * MM_PER_INCH
    if bounds.zmin < -thickness - 1e-5 or bounds.zmax > 1e-5:
        raise RuntimeError("Machining changed the plate's vertical bounds")
    final_volume = solid_volume(plate)
    if final_volume > original_volume + 1e-5 or final_volume <= 0:
        raise RuntimeError("Invalid material-removal volume")
    if active and original_volume - final_volume <= 1e-7:
        raise RuntimeError("Cutter model did not remove any material")
    return plate, {
        **model_report,
        "active_grooves": len(active),
        "solid_valid": True,
        "solid_count": 1,
        "faces": len(plate.Faces()),
        "original_volume_mm3": original_volume,
        "final_volume_mm3": final_volume,
        "removed_volume_mm3": original_volume - final_volume,
        "minimum_remaining_thickness_mm": thickness - parameters.max_depth * MM_PER_INCH,
    }


def export_obj(shape, path: Path, tolerance_mm=0.05):
    if not math.isfinite(tolerance_mm) or tolerance_mm <= 0:
        raise ValueError("OBJ mesh tolerance must be finite and positive (mm)")
    # CadQuery 2.5's default mesher uses relative deflection. Mesh explicitly
    # with absolute millimeter deflection, then extract its triangulation.
    BRepTools.Clean_s(shape.wrapped)
    mesher = BRepMesh_IncrementalMesh(shape.wrapped, tolerance_mm, False, 0.1, True)
    if not mesher.IsDone():
        raise RuntimeError("Finished plate could not be meshed")
    vertices, triangles = shape.tessellate(tolerance_mm, 0.1)
    if not vertices or not triangles:
        raise RuntimeError("Finished plate could not be tessellated")
    # Weld coincident CAD face-boundary vertices. OBJ coordinates retain the
    # same millimeter convention as STEP (OBJ itself has no unit metadata).
    unique, indices = [], {}
    mapping = []
    for vertex in vertices:
        key = tuple(round(value, 7) for value in vertex.toTuple())
        if key not in indices:
            indices[key] = len(unique) + 1
            unique.append(vertex.toTuple())
        mapping.append(indices[key])
    faces = []
    for triangle in triangles:
        face = tuple(mapping[index] for index in triangle)
        if len(set(face)) == 3:
            faces.append(face)
    edges = Counter((a, b) for face in faces
                    for a, b in zip(face, (face[1], face[2], face[0])))
    if not faces or any(count != 1 or edges.get((b, a), 0) != 1
                        for (a, b), count in edges.items()):
        raise RuntimeError("OBJ tessellation is not a closed, consistently oriented mesh")
    points = np.asarray(unique)
    face_indices = np.asarray(faces) - 1
    a, b, c = (points[face_indices[:, i]] for i in range(3))
    mesh_volume = float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6)
    if mesh_volume <= 0:
        raise RuntimeError("OBJ tessellation has invalid orientation or volume")
    with path.open("w", encoding="utf-8") as output:
        output.write("# MillCraft machined plate; coordinates in millimeters\no machined_plate\n")
        for x, y, z in unique:
            output.write(f"v {x:.12g} {y:.12g} {z:.12g}\n")
        for a, b, c in faces:
            output.write(f"f {a} {b} {c}\n")
    return {"vertices": len(unique), "triangles": len(faces), "tolerance_mm": tolerance_mm,
            "watertight": True, "signed_volume_mm3": mesh_volume}
