"""Portable SVG plot of centerline depth, not a cutter/material simulation."""
from html import escape
from pathlib import Path

import numpy as np

from .sampling import MM_PER_INCH


def write_preview(rows, parameters, path: Path):
    w = parameters.artwork_width * MM_PER_INCH
    h = parameters.artwork_height * MM_PER_INCH
    depth = parameters.max_depth * MM_PER_INCH
    minimum_depth = parameters.min_depth * MM_PER_INCH
    depth_range = depth - minimum_depth
    plot_w, plot_h = 520 * w / max(w, h), 520 * h / max(w, h)
    left, top = 50 + (520 - plot_w) / 2, 120 + (520 - plot_h) / 2
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="740" viewBox="0 0 1100 740">',
             '<rect width="1100" height="740" fill="#f7f7f4"/>',
             '<g font-family="Arial, sans-serif" fill="#222">']

    def label(x, y, text, size=16):
        parts.append(f'<text x="{x}" y="{y}" font-size="{size}">{escape(text)}</text>')

    label(50, 46, "MillCraft / spline depth preview", 26)
    label(50, 77, f"Centerline geometry only · white = {parameters.min_depth:g} in depth · black = {parameters.max_depth:g} in depth")
    parts.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="white" stroke="#b5b5b0"/>')
    stroke = max(0.5, min(8, parameters.line_spacing * MM_PER_INCH / h * plot_h * 0.65))
    for row in rows:
        x = np.linspace(-w / 2, w / 2, 601)
        z = row.polynomial((x[:-1] + x[1:]) / 2)
        y = top + (0.5 - row.y_mm / h) * plot_h
        for a, b, zz in zip(x[:-1], x[1:], z):
            relative_depth = (-zz - minimum_depth) / depth_range if depth_range > 0 else 0
            shade = int(round(245 * (1 - np.clip(relative_depth, 0, 1))))
            color = f"rgb({shade},{shade},{shade})"
            x1, x2 = left + (a / w + 0.5) * plot_w, left + (b / w + 0.5) * plot_w
            parts.append(f'<path d="M{x1:.3f},{y:.3f}H{x2:.3f}" stroke="{color}" stroke-width="{stroke:.3f}"/>')
    label(50, 675, f"Top view · {parameters.artwork_width:g} × {parameters.artwork_height:g} in · +Y up")
    label(620, 140, "Representative row profiles", 20)
    graph_x, graph_y, graph_w, graph_h = 635, 185, 390, 210
    for fraction in (0, 0.5, 1):
        yy = graph_y + graph_h * fraction
        parts.append(f'<path d="M{graph_x},{yy}h{graph_w}" stroke="#d0d0ca"/>')
        label(graph_x, yy - 8, f"Z = {-parameters.max_depth * fraction:.4f} in", 12)
    indices = sorted(set([len(rows) // 4, len(rows) // 2, 3 * len(rows) // 4]))
    for j, (index, color) in enumerate(zip(indices, ("#b3611f", "#185f99", "#45802d"))):
        row = rows[index]
        x = np.linspace(-w / 2, w / 2, 601)
        z = row.polynomial(x)
        points = " ".join(f"{graph_x + (xx / w + 0.5) * graph_w:.3f},{graph_y - zz / depth * graph_h:.3f}"
                          for xx, zz in zip(x, z))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>')
        parts.append(f'<rect x="635" y="{435 + 28 * j}" width="18" height="3" fill="{color}"/>')
        label(665, 441 + 28 * j, f"Y = {row.y_mm / MM_PER_INCH:+.3f} in", 14)
    label(635, 420, f"X = {-parameters.artwork_width / 2:g} to {parameters.artwork_width / 2:g} in", 13)
    label(620, 562, f"{len(rows)} rows · {len(rows[0].polynomial.x)} samples / row")
    label(620, 592, f"Row spacing: {parameters.line_spacing:g} in")
    label(620, 622, f"Sample spacing: {parameters.sample_spacing:g} in")
    label(620, 652, f"Depth range: {parameters.min_depth:g}–{parameters.max_depth:g} in")
    label(50, 712, "CAD coordinates: mm · plate centered at XY=(0,0) · plate top Z=0 · depths negative", 14)
    parts.append('</g></svg>')
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")
