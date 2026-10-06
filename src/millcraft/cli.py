import argparse
from dataclasses import asdict
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

from .preview import write_preview
from .sampling import MM_PER_INCH, Parameters, load_grayscale, sample_image
from .splines import interpolate_grid


def parser():
    result = argparse.ArgumentParser(
        description="Cut image-derived grooves into a plate and export STEP + OBJ (mm).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument("image", type=Path, help="Input image; use image.png for the supplied blurred A")
    result.add_argument("--out-dir", type=Path, default=Path("outputs/a"))
    defaults = Parameters()
    for name, value in asdict(defaults).items():
        result.add_argument("--" + name.replace("_", "-"), type=float, default=value,
                            help=f"{name.replace('_', ' ').capitalize()} in inches")
    result.add_argument("--image-fit", choices=("contain", "stretch"), default="contain",
                        help="Preserve image aspect ratio with white margins, or stretch to artwork")
    result.add_argument("--blur-px", type=float, default=0,
                        help="Optional additional Gaussian blur in source-image pixels")
    result.add_argument("--spline", choices=("pchip", "cubic"), default="pchip",
                        help="C1 shape-preserving or C2 natural cubic; reject Z overshoot in either mode")
    result.add_argument("--with-plate", action="store_true", help="In centerlines mode, also export a separate uncut plate")
    result.add_argument("--geometry", choices=("machined", "centerlines"), default="machined",
                        help="Finished plate solid, or the original curve-only workflow")
    result.add_argument("--cut-sample-spacing", type=float, default=0.01,
                        help="Maximum X spacing for the cutter envelope in inches")
    result.add_argument("--cutter-segments", type=int, default=32,
                        help="Circular cutter polygon resolution; multiple of 4, from 8 to 128")
    result.add_argument("--mesh-tolerance-mm", type=float, default=0.05,
                        help="OBJ tessellation deflection in millimeters")
    return result


def generate(args, progress=None):
    parameters = Parameters(**{name: getattr(args, name) for name in asdict(Parameters())})
    gray = load_grayscale(args.image, args.blur_px)
    grid = sample_image(gray, parameters, args.image_fit)
    rows = interpolate_grid(grid, parameters.max_depth * MM_PER_INCH, args.spline,
                            min_depth_mm=parameters.min_depth * MM_PER_INCH)
    # Import CAD only after image sampling and analytic validation succeed.
    from .step import export_centerlines, export_plate

    machined = args.geometry == "machined"
    outputs = (["machined_plate.step", "machined_plate.obj"] if machined else ["centerlines.step"])
    outputs += ["samples.csv", "report.json", "preview.svg"]
    if args.with_plate and not machined:
        outputs.append("plate_reference.step")
    source_path = args.image.resolve()
    if any((args.out_dir / name).resolve() == source_path for name in outputs):
        raise ValueError("An output path would overwrite the source image")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    obsolete = []
    # Retire only known generated files listed by the preceding run. An
    # unrelated STEP file placed in this folder by the user is preserved.
    previous_report = args.out_dir / "report.json"
    if previous_report.is_file():
        try:
            previous_outputs = json.loads(previous_report.read_text(encoding="utf-8")).get("outputs", [])
            known_geometry = {"centerlines.step", "plate_reference.step", "machined_plate.step", "machined_plate.obj"}
            obsolete = [name for name in previous_outputs
                        if name in known_geometry and name not in outputs
                        and (args.out_dir / name).resolve() != source_path]
        except (ValueError, OSError, AttributeError, TypeError):
            pass
    # Stage all artifacts so a failed validation/export does not replace a good STEP.
    with TemporaryDirectory(prefix=".millcraft-", dir=args.out_dir) as staging:
        stage = Path(staging)
        machining_report, mesh_report = None, None
        if machined:
            import cadquery as cq
            from .machining import export_obj, machine_plate

            solid, machining_report = machine_plate(
                rows, parameters, args.cut_sample_spacing, args.cutter_segments, progress)
            if progress:
                progress("Exporting finished plate STEP and OBJ...")
            cq.exporters.export(solid, str(stage / "machined_plate.step"), exportType="STEP")
            mesh_report = export_obj(solid, stage / "machined_plate.obj", args.mesh_tolerance_mm)
        else:
            export_centerlines(rows, stage / "centerlines.step")
        if args.with_plate and not machined:
            export_plate(parameters, stage / "plate_reference.step")
        xx, yy = np.meshgrid(grid.x_mm, grid.y_mm)
        points = np.column_stack((np.repeat(np.arange(len(rows)), len(grid.x_mm)),
                                  xx.ravel(), yy.ravel(), grid.z_mm.ravel()))
        np.savetxt(stage / "samples.csv", points, delimiter=",", header="row,x_mm,y_mm,z_mm",
                   comments="", fmt=["%d", "%.12f", "%.12f", "%.12f"])
        report = {
            "geometry": "machined_plate_solid" if machined else "spline_centerlines_only",
            "source": {"path": str(source_path), "sha256": hashlib.sha256(args.image.read_bytes()).hexdigest(),
                       "size_px": list(grid.image_size)},
            "units": {"input_parameters": "inch", "geometry": "mm", "samples_csv": "mm", "step": "mm", "obj": "mm"},
            "coordinates": "plate centered in XY; top Z=0; image top is +Y; darker is more negative Z",
            "plate_reference": {"material": "aluminum", "machined": machined},
            "machining": machining_report,
            "obj_mesh": mesh_report,
            "parameters_inch": asdict(parameters),
            "parameters_mm": {k: v * MM_PER_INCH for k, v in asdict(parameters).items()},
            "processing": {"image_fit": args.image_fit, "additional_blur_px": args.blur_px,
                           "depth_mapping": "z_mm = -(min_depth_mm + (1 - grayscale / 255) * (max_depth_mm - min_depth_mm))",
                           "spline": args.spline, "continuity": "C1" if args.spline == "pchip" else "C2"},
            "grid": {"rows": len(rows), "samples_per_row": len(grid.x_mm), "total_samples": grid.z_mm.size,
                     "x_extent_mm": [float(grid.x_mm[0]), float(grid.x_mm[-1])],
                     "y_extent_mm": [float(grid.y_mm[0]), float(grid.y_mm[-1])],
                     "final_x_spacing_mm": float(grid.x_mm[-1] - grid.x_mm[-2]),
                     "final_y_spacing_mm": float(grid.y_mm[-1] - grid.y_mm[-2])},
            "validation": {"passed": True, "method": "endpoints and all real cubic derivative roots per span",
                           "tolerance_mm": 1e-7,
                           "z_min_mm": min(row.z_min_mm for row in rows),
                           "z_max_mm": max(row.z_max_mm for row in rows),
                           "maximum_local_overshoot_mm": max(row.local_overshoot_mm for row in rows),
                           "global_depth_overshoot_mm": max(0.0, -parameters.max_depth * MM_PER_INCH
                                                            - min(row.z_min_mm for row in rows)),
                           "above_surface_overshoot_mm": max(0.0, max(row.z_max_mm for row in rows)),
                           "minimum_depth_undershoot_mm": max(0.0, parameters.min_depth * MM_PER_INCH
                                                               + max(row.z_max_mm for row in rows)),
                           "occt_conversion_checked": not machined},
            "rows": [{"index": i, "y_mm": row.y_mm, "z_min_mm": row.z_min_mm,
                      "z_max_mm": row.z_max_mm, "local_overshoot_mm": row.local_overshoot_mm}
                     for i, row in enumerate(rows)],
            "tool": {"diameter_inch": parameters.tool_diameter,
                     "use": ("vertical flat end mill swept-envelope material removal; no G-code"
                             if machined else "centerline metadata only")},
            "outputs": outputs,
            "packages": {name: version(name) for name in ("millcraft", "Pillow", "numpy", "scipy", "cadquery", "cadquery-ocp")},
        }
        (stage / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        write_preview(rows, parameters, stage / "preview.svg")
        for name in outputs:
            (stage / name).replace(args.out_dir / name)
        for name in obsolete:
            (args.out_dir / name).unlink(missing_ok=True)
    return report


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        report = generate(args, progress=lambda message: print(message, flush=True))
    except (ValueError, OSError, RuntimeError) as error:
        print(f"millcraft: {error}", file=sys.stderr)
        return 1
    grid, validation = report["grid"], report["validation"]
    if args.geometry == "machined":
        print(f"Finished plate: {args.out_dir / 'machined_plate.step'}")
        print(f"OBJ preview: {args.out_dir / 'machined_plate.obj'}")
        print(f"Removed volume: {report['machining']['removed_volume_mm3']:.3f} mm³; one valid solid")
    else:
        print(f"Exported {grid['rows']} spline rows × {grid['samples_per_row']} samples to {args.out_dir / 'centerlines.step'}")
    print(f"Analytic Z range: [{validation['z_min_mm']:.9f}, {validation['z_max_mm']:.9f}] mm; validation passed")
    if args.with_plate and args.geometry == "centerlines":
        print(f"Plate reference: {args.out_dir / 'plate_reference.step'}")
    print(f"Preview: {args.out_dir / 'preview.svg'}; samples and validation report included")
    return 0
