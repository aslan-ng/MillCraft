"""Translate validated curves to exact OpenCascade B-splines and STEP."""
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCP.Geom import Geom_BSplineCurve
from OCP.TColgp import TColgp_Array1OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.gp import gp_Pnt

from .sampling import MM_PER_INCH


def row_edge(row):
    """One cubic BSpline edge per horizontal row, exactly matching SciPy.

    Convert each cubic Hermite span to Bezier poles, then concatenate spans
    into a nonrational B-spline (degree 3, internal knot multiplicity 3).
    PCHIP remains C1; natural cubic remains C2 despite knot multiplicity.
    """
    curve = row.polynomial
    x, z = curve.x, curve(curve.x)
    dz = curve.derivative()(x)
    poles = [(x[0], row.y_mm, z[0])]
    for i, h in enumerate(np.diff(x)):
        poles.extend([
            (x[i] + h / 3, row.y_mm, z[i] + h * dz[i] / 3),
            (x[i + 1] - h / 3, row.y_mm, z[i + 1] - h * dz[i + 1] / 3),
            (x[i + 1], row.y_mm, z[i + 1]),
        ])
    points = TColgp_Array1OfPnt(1, len(poles))
    knots = TColStd_Array1OfReal(1, len(x))
    multiplicities = TColStd_Array1OfInteger(1, len(x))
    for i, p in enumerate(poles, 1):
        points.SetValue(i, gp_Pnt(*map(float, p)))
    for i, knot in enumerate(x, 1):
        knots.SetValue(i, float(knot))
        multiplicities.SetValue(i, 4 if i in (1, len(x)) else 3)
    spline = Geom_BSplineCurve(points, knots, multiplicities, 3, False)
    builder = BRepBuilderAPI_MakeEdge(spline)
    if not builder.IsDone():
        raise RuntimeError(f"OpenCascade could not construct row Y={row.y_mm} mm")
    edge = cq.Edge(builder.Edge())
    if not edge.isValid():
        raise RuntimeError(f"Invalid spline edge at Y={row.y_mm} mm")
    # Verify conversion at all analytic extrema and extra interior points.
    interior = (x[:-1, None] + np.diff(x)[:, None] * np.array([0.25, 0.5, 0.75])).ravel()
    for t in np.unique(np.concatenate((x, row.extrema_x_mm, interior))):
        p = spline.Value(float(t))
        if max(abs(p.X() - t), abs(p.Y() - row.y_mm),
               abs(p.Z() - float(curve(t)))) > 1e-7:
            raise RuntimeError("OpenCascade conversion changed the interpolated geometry")
    return edge


def export_centerlines(rows, path: Path):
    edges = [row_edge(row) for row in rows]
    if not edges:
        raise ValueError("At least one row is required")
    cq.exporters.export(cq.Compound.makeCompound(edges), str(path), exportType="STEP")
    return edges


def plate_shape(parameters):
    # Centered in XY, top at Z=0, underside at -thickness.
    return cq.Workplane("XY").box(
        parameters.plate_width * MM_PER_INCH,
        parameters.plate_height * MM_PER_INCH,
        parameters.plate_thickness * MM_PER_INCH,
        centered=(True, True, False),
    ).translate((0, 0, -parameters.plate_thickness * MM_PER_INCH)).val()


def export_plate(parameters, path: Path):
    cq.exporters.export(plate_shape(parameters), str(path), exportType="STEP")
