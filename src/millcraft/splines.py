"""Exact piecewise cubic interpolation and analytical Z-extrema checks."""
from dataclasses import dataclass
import math

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, PPoly


class SplineOvershootError(ValueError):
    pass


@dataclass
class RowSpline:
    y_mm: float
    polynomial: PPoly
    z_min_mm: float
    z_max_mm: float
    local_overshoot_mm: float
    extrema_x_mm: np.ndarray


def interval_extrema(coefficients: np.ndarray, length: float) -> np.ndarray:
    """Endpoints plus every real derivative root inside a cubic interval.

    Solve in normalized t=[0,1] so short/long sample spacing does not
    distort the derivative's coefficient magnitudes.
    """
    a, b, c, _ = coefficients
    roots = np.roots([3 * a * length**3, 2 * b * length**2, c * length])
    interior = [float(r.real) for r in roots
                if abs(r.imag) < 1e-10 and 0 < r.real < 1]
    return np.array([0.0, *interior, 1.0]) * length


def interpolate_row(x_mm, y_mm, z_mm, max_depth_mm, method="pchip", tolerance_mm=1e-7,
                    min_depth_mm=0.0):
    x, z = np.asarray(x_mm, dtype=float), np.asarray(z_mm, dtype=float)
    if (x.ndim != 1 or z.shape != x.shape or len(x) < 2
            or not np.isfinite(x).all() or not np.isfinite(z).all()
            or np.any(np.diff(x) <= 0)):
        raise ValueError("Row needs finite Z values and strictly increasing X samples")
    if not math.isfinite(y_mm) or not math.isfinite(max_depth_mm) or max_depth_mm <= 0:
        raise ValueError("Row Y and positive maximum depth must be finite")
    if not math.isfinite(tolerance_mm) or tolerance_mm < 0:
        raise ValueError("Spline tolerance must be finite and nonnegative")
    if not math.isfinite(min_depth_mm) or not 0 <= min_depth_mm <= max_depth_mm:
        raise ValueError("Minimum depth must be finite and between zero and maximum depth")
    if method == "pchip":
        curve = PchipInterpolator(x, z, extrapolate=False)
    elif method == "cubic":
        curve = CubicSpline(x, z, bc_type="natural", extrapolate=False)
    else:
        raise ValueError("Spline method must be 'pchip' or 'cubic'")
    extrema, values, local_overshoot = [], [], 0.0
    for i, h in enumerate(np.diff(x)):
        dx = interval_extrema(curve.c[:, i], h)
        local_values = np.polyval(curve.c[:, i], dx)
        extrema.extend(x[i] + dx)
        values.extend(local_values)
        lo, hi = min(z[i:i + 2]), max(z[i:i + 2])
        local_overshoot = max(local_overshoot, lo - min(local_values), max(local_values) - hi)
    z_min, z_max = float(min(values)), float(max(values))
    if z_min < -max_depth_mm - tolerance_mm or z_max > -min_depth_mm + tolerance_mm:
        raise SplineOvershootError(
            f"Row Y={y_mm:.6f} mm: spline Z=[{z_min:.9f}, {z_max:.9f}] mm "
            f"exceeds [-{max_depth_mm:.9f}, {-min_depth_mm:.9f}] mm. Use --spline pchip."
        )
    if local_overshoot > tolerance_mm:
        raise SplineOvershootError(
            f"Row Y={y_mm:.6f} mm: spline overshoots adjacent sample depths by "
            f"{local_overshoot:.9f} mm. Use --spline pchip."
        )
    return RowSpline(float(y_mm), curve, z_min, z_max,
                     float(local_overshoot), np.unique(extrema))


def interpolate_grid(grid, max_depth_mm, method="pchip", min_depth_mm=0.0):
    return [interpolate_row(grid.x_mm, y, z, max_depth_mm, method, min_depth_mm=min_depth_mm)
            for y, z in zip(grid.y_mm, grid.z_mm)]
