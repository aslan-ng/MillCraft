"""Run MillCraft from VS Code; edit the settings below before running."""

from pathlib import Path

from millcraft.cli import main as cli_main


PROJECT_DIR = Path(__file__).resolve().parent

# Use a project-relative path or replace it with an absolute path.
IMAGE = PROJECT_DIR / "image.png"
OUTPUT_DIR = PROJECT_DIR / "outputs" / "custom"

# All dimensions below are in inches. Exported geometry uses millimeters.
PARAMETERS = {
    "plate_width": 3.5,
    "plate_height": 3.5,
    "plate_thickness": 0.25,
    "artwork_width": 3.0,
    "artwork_height": 3.0,
    "tool_diameter": 1 / 8 + 0.01,
    "line_spacing": 0.15,
    "sample_spacing": 0.025,
    "min_depth": 0.03,  # Pure white depth; positive values also cut white areas.
    "max_depth": 0.19,
    "image_fit": "contain",  # "contain" or "stretch"
    "blur_px": 0,  # Additional blur radius in source-image pixels.
    "spline": "pchip",  # "pchip" or "cubic"; cubic can fail overshoot validation.
}
GEOMETRY = "machined"  # Finished plate STEP + OBJ; "centerlines" exports curves.
CUT_MODEL = "spline-pocket"  # Rectangular rows with one spline floor per row.
MESH_TOLERANCE_MM = 0.05  # OBJ detail; smaller gives a finer mesh.


def main():
    args = [str(IMAGE), "--out-dir", str(OUTPUT_DIR), "--geometry", GEOMETRY,
            "--cut-model", CUT_MODEL,
            "--mesh-tolerance-mm", str(MESH_TOLERANCE_MM)]
    for name, value in PARAMETERS.items():
        args.extend(["--" + name.replace("_", "-"), str(value)])
    return cli_main(args)


if __name__ == "__main__":
    raise SystemExit(main())
