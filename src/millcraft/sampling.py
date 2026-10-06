"""Image sampling and physical coordinates; no CAD or cutter operations."""
from dataclasses import asdict, dataclass
from pathlib import Path
import math

import numpy as np
from PIL import Image, ImageFilter, ImageOps

MM_PER_INCH = 25.4


@dataclass(frozen=True)
class Parameters:
    # Public inputs are inches. Arrays and CAD geometry are millimeters.
    plate_width: float = 3.5
    plate_height: float = 3.5
    plate_thickness: float = 0.25
    artwork_width: float = 3.0
    artwork_height: float = 3.0
    tool_diameter: float = 1 / 16
    line_spacing: float = 0.075
    sample_spacing: float = 0.025
    max_depth: float = 0.035
    min_depth: float = 0.0

    def __post_init__(self):
        for name, value in asdict(self).items():
            if name == "min_depth":
                if not math.isfinite(value) or value < 0:
                    raise ValueError("Minimum depth must be finite and nonnegative (inches)")
                continue
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive (inches)")
        if self.artwork_width > self.plate_width or self.artwork_height > self.plate_height:
            raise ValueError("Artwork must fit inside the plate")
        if self.max_depth >= self.plate_thickness:
            raise ValueError("Maximum depth must be less than plate thickness")
        if self.min_depth > self.max_depth:
            raise ValueError("Minimum depth must not exceed maximum depth")
        if ((math.ceil(self.artwork_width / self.sample_spacing) + 1)
                * (math.ceil(self.artwork_height / self.line_spacing) + 1) > 2_000_000):
            raise ValueError("Requested grid exceeds the 2,000,000 sample limit")


@dataclass
class SampleGrid:
    x_mm: np.ndarray
    y_mm: np.ndarray
    z_mm: np.ndarray  # [row, column], rows run bottom to top
    image_size: tuple


def inclusive_axis(length_mm: float, spacing_mm: float) -> np.ndarray:
    """Keep nominal spacing, with a shorter final interval when necessary."""
    if not all(math.isfinite(v) and v > 0 for v in (length_mm, spacing_mm)):
        raise ValueError("Axis length and spacing must be finite and positive")
    count = math.floor(length_mm / spacing_mm)
    offsets = np.arange(count + 1, dtype=float) * spacing_mm
    if math.isclose(float(offsets[-1]), length_mm, abs_tol=1e-9, rel_tol=1e-12):
        offsets[-1] = length_mm
    else:
        offsets = np.append(offsets, length_mm)
    return offsets - length_mm / 2


def load_grayscale(path: Path, blur_px: float = 0) -> np.ndarray:
    if not math.isfinite(blur_px) or blur_px < 0:
        raise ValueError("Blur radius must be finite and nonnegative")
    with Image.open(path) as source:
        oriented = ImageOps.exif_transpose(source).convert("RGBA")
        white = Image.new("RGBA", oriented.size, "white")
        white.alpha_composite(oriented)
        gray = white.convert("L")
        if blur_px:
            gray = gray.filter(ImageFilter.GaussianBlur(blur_px))
        return np.asarray(gray, dtype=float) / 255.0


def sample_image(gray: np.ndarray, parameters: Parameters, fit: str = "contain") -> SampleGrid:
    """Bilinear samples: white -> Z=-min_depth, black -> Z=-max_depth.

    Image top maps to positive Y; X/Y are centered on the plate. `contain`
    preserves image aspect ratio and uses white outside its footprint.
    """
    gray = np.asarray(gray, dtype=float)
    if gray.ndim != 2 or min(gray.shape) < 2:
        raise ValueError("Image must be a grayscale array at least 2 by 2 pixels")
    if not np.isfinite(gray).all() or np.any((gray < 0) | (gray > 1)):
        raise ValueError("Grayscale values must be finite and in [0, 1]")
    if fit not in ("contain", "stretch"):
        raise ValueError("Image fit must be 'contain' or 'stretch'")
    width = parameters.artwork_width * MM_PER_INCH
    height = parameters.artwork_height * MM_PER_INCH
    x = inclusive_axis(width, parameters.sample_spacing * MM_PER_INCH)
    y = inclusive_axis(height, parameters.line_spacing * MM_PER_INCH)
    image_h, image_w = gray.shape
    if fit == "contain":
        scale = min(width / image_w, height / image_h)
        footprint_w, footprint_h = image_w * scale, image_h * scale
    else:
        footprint_w, footprint_h = width, height
    xx, yy = np.meshgrid(x, y)
    u = (xx / footprint_w + 0.5) * (image_w - 1)
    v = (0.5 - yy / footprint_h) * (image_h - 1)
    inside = ((u >= -1e-10) & (u <= image_w - 1 + 1e-10)
              & (v >= -1e-10) & (v <= image_h - 1 + 1e-10))
    u, v = np.clip(u, 0, image_w - 1), np.clip(v, 0, image_h - 1)
    left, top = np.floor(u).astype(int), np.floor(v).astype(int)
    right, bottom = np.minimum(left + 1, image_w - 1), np.minimum(top + 1, image_h - 1)
    dx, dy = u - left, v - top
    intensity = ((1 - dx) * (1 - dy) * gray[top, left]
                 + dx * (1 - dy) * gray[top, right]
                 + (1 - dx) * dy * gray[bottom, left]
                 + dx * dy * gray[bottom, right])
    intensity = np.where(inside, intensity, 1.0)
    z = -(parameters.min_depth + (1.0 - intensity)
          * (parameters.max_depth - parameters.min_depth)) * MM_PER_INCH
    return SampleGrid(x, y, z, (image_w, image_h))
