# MillCraft

Turn a blurred grayscale image into horizontal toolpaths, cut their flat end
mill grooves into a plate, and export the finished solid as STEP and OBJ.
The supplied `image.png` is a 2100 × 2100 blurred A;
it is used directly, with no additional blur by default. The original attachment
name `A-flat-blurred.png` is not present in this checkout.

The default output is one `machined_plate.step` solid and a matching
`machined_plate.obj` mesh for inspection. The original curve-only export is
available with `--geometry centerlines`.

## Run locally

Use Python 3.9–3.12. Dependencies are open-source Pillow, NumPy, SciPy, and
CadQuery/OpenCascade. CadQuery 2.5.2 is pinned for compatibility with the Python
3.9 / Apple Silicon environment used for validation.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m millcraft image.png --out-dir outputs/a
```

The virtual environment is already installed in this workspace. To regenerate
the current artifacts, run only the last command. The installed `millcraft`
command is also available at `.venv/bin/millcraft`.

### Run as a Python file in VS Code

Open `main.py` in the project root. Edit `IMAGE`, `OUTPUT_DIR`, `PARAMETERS`,
and the machining settings to choose your input and settings. All dimensions in
`PARAMETERS` are inches; `blur_px` is in source-image pixels.

In VS Code, use **Python: Select Interpreter** from the Command Palette and
select this project's `.venv/bin/python`. Then use **Run Python File in Terminal**
with `main.py` open, or run it from the integrated terminal:

```sh
.venv/bin/python main.py
```

The default script settings write to `outputs/custom/`. Paths are relative to
the script's directory, so the script also works from another working directory.
Rerunning replaces the generated files in the chosen output directory.
Artwork must fit inside the plate, and maximum depth must be less than plate
thickness.

| Parameter | Default, inches | Millimeters |
| --- | ---: | ---: |
| Aluminum plate width × height | 3.5 × 3.5 | 88.9 × 88.9 |
| Plate thickness | 0.25 | 6.35 |
| Artwork width × height | 3.0 × 3.0 | 76.2 × 76.2 |
| Primary flat end mill diameter | 0.0625 (1/16) | 1.5875 |
| Horizontal row spacing | 0.075 | 1.905 |
| Samples along each row | 0.025 spacing | 0.635 spacing |
| Minimum depth (pure white) | 0.0 | 0.0 |
| Maximum depth | 0.035 | 0.889 |

All dimensional CLI arguments use inches. CAD geometry, STEP units, and sample
coordinates use millimeters. The defaults produce **41 rows × 121 samples =
4,961 points** and **41 spline edges**. Both artwork boundaries are included.
When a dimension is not divisible by its spacing, the final interval is shorter;
the preceding intervals keep the specified spacing.

For another input or parameter set:

```sh
.venv/bin/python -m millcraft A-flat-blurred.png \
  --out-dir outputs/custom \
  --line-spacing 0.075 --sample-spacing 0.025 --max-depth 0.035 \
  --geometry machined
.venv/bin/python -m millcraft --help
```

The artwork is centered on the plate, giving a 0.25 inch margin on each side.
XY=(0,0) is the plate center. The top face is Z=0, and the underside is Z=-6.35 mm.
Image left is negative X, and image top is positive Y. Rows run left to right,
ordered from the bottom to the top of the image.

## Image and depth mapping

1. Apply EXIF orientation, composite transparency onto white, and convert to
   grayscale with Pillow.
2. Optionally apply `--blur-px RADIUS` in source-image pixels. Default: 0, because
   the supplied A is already blurred.
3. Map the image to the artwork region. `--image-fit contain` preserves its
   aspect ratio, centers it, and fills unused space with white. `stretch` fills
   the artwork by scaling each axis independently. Neither mode crops the image.
4. Bilinearly sample grayscale at the physical grid coordinates. Pixel indices
   0 and size-1 map to the edges of the image footprint.
5. Map normalized intensity `g` linearly to depth:

   ```text
   Z_mm = -(minimum_depth_mm + (1 - g) × (maximum_depth_mm - minimum_depth_mm))
   white (255) -> 0
   gray  (127.5, interpolated) -> -0.4445 mm
   black (0) -> -0.889 mm
   ```

   Set `"min_depth"` in `main.py` or pass `--min-depth` to change the pure-white
   depth. Inputs are positive depths in inches; CAD Z is negative. For example,
   `min_depth=0.01` and `max_depth=0.20` map white to a 0.01-inch cut,
   50% gray to 0.105 inches, and black to 0.20 inches. Minimum depth must be
   finite and satisfy `0 <= min_depth <= max_depth`. Equal limits make a
   constant-depth path. White containment margins and transparency also map
   to minimum depth. With positive minimum depth, white rows are cut as well.

No contrast normalization or extra smoothing is applied by default. White rows
are retained as centerlines at Z=-min_depth (zero by default).

## Splines and overshoot validation

The default `--spline pchip` uses SciPy's shape-preserving piecewise cubic Hermite
interpolator. It passes through every sampled point with **C1 continuity**
(continuous position and tangent); its second derivative can change at sample
knots. Shape preservation avoids introducing new depth peaks between samples.
See the [SciPy PCHIP documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html).

`--spline cubic` requests a natural **C2** cubic interpolator. This can overshoot,
especially near a flat black or white region, so export fails if validation
detects an overshoot. It does not silently clamp depths or change spline methods.

For **every span**, validation evaluates its endpoints and all real roots of
the cubic's Z derivative inside the span. This finds the actual Z extrema rather
than relying on a dense sampling approximation. Both the global depth range
`[-max_depth, -min_depth]` and the adjacent samples' depth range must be respected within
`1e-7 mm`. Tiny floating-point residues remain visible in the report.

In `--geometry centerlines` mode, the validated cubics are translated exactly into Bezier control points, then
assembled into **one degree-3 OpenCascade B-spline edge per row**. Internal knots
have multiplicity 3; the control points preserve the interpolator's actual C1
or C2 continuity. X stays monotonic and Y stays constant. Before export, the
OpenCascade curve is checked against the interpolator at all sample knots,
all Z extrema, and three interior points in each span. There is no additional
fitting step that could change the depth curve.

STEP export uses CadQuery's exporter and its millimeter convention. See the
[CadQuery API documentation](https://cadquery.readthedocs.io/en/stable/classreference.html).

## Cut model and accuracy

The default `--cut-model spline-pocket` makes rectangular row pockets. Each row
has one exact spline defining its floor, extruded uniformly across the width.
At any X position, the floor has the same Z everywhere across that row's width.
The top outline has straight sides and straight ends at the spline's X endpoints.
Pocket width is `tool_diameter`; there are no rounded end caps or extensions
beyond the row endpoints. In `main.py`, `CUT_MODEL = "spline-pocket"` selects
this model.

Each pocket volume has one curved floor, four planar walls, and a planar top.
The STEP stores that floor as a single extrusion surface rather than many
transverse cutter sections. This gives much smaller CAD files. The floor
matches the analytically validated row spline, including its minimum/maximum
depth settings, without another fitting or resampling step. Rows whose depth is
zero remain uncut. A positive `min_depth` cuts all white portions as well.
Overlapping pockets remove their combined volumes.

`MESH_TOLERANCE_MM` / `--mesh-tolerance-mm` controls OBJ tessellation deflection
in millimeters (default 0.05). Smaller values make a finer OBJ mesh; the STEP
floor stays the same. `report.json` records pocket dimensions, face count,
material removed, and the closed mesh's triangle count. Volume reports use
adaptive integration over the spline spans.

The output is checked for one valid solid and positive remaining volume. The
floor stays within the requested depth range, preserving at least plate thickness
minus maximum depth. The plate edges clip any pocket extending outside the stock.
This is a direct model of the requested pocket geometry.

### Optional circular cutter model

The previous rounded cutter model remains available with
`--cut-model cutter-envelope`. It approximates the swept envelope of a vertical,
cylindrical flat end mill. At each transverse section, the floor is the lowest
spline depth under the cutter's circular footprint. It includes the finite-radius
effect on slopes and rounded ends, using linear X samples and ruled sections
around an inscribed polygonal cutter.

Only in this mode, `--cut-sample-spacing` sets envelope X spacing in inches
(default 0.01), and `--cutter-segments` sets circular footprint resolution
(default 32; a multiple of 4 between 8 and 128). The reported radial chord error
is not a bound on total floor error. Near Z=0, a maximum 0.0001 mm upward clearance
avoids grazing CAD faces; rows shallower than that clearance are skipped. Floor
samples deeper than 0.001 mm have no clearance offset. This mode is an
approximation and generally creates larger files.

Neither mode generates lead-ins, linking moves, feeds/speeds, or G-code.

## Output files

| File, under `outputs/a/` | Contents |
| --- | --- |
| `machined_plate.step` | One solid plate with the grooves cut into it |
| `machined_plate.obj` | Triangle mesh tessellated from the same finished solid |
| `samples.csv` | Sample coordinates: `row,x_mm,y_mm,z_mm` |
| `report.json` | Parameters, units, spline validation, pocket dimensions, removed volume, and mesh counts |
| `preview.svg` | Top view colored by spline depth and representative Z profiles |

The default run produces one STEP file. STEP and OBJ use the same coordinates
in millimeters; OBJ has no standard unit metadata, so choose millimeters if an
importer asks. Open the OBJ in a mesh viewer before importing the STEP into CAD.
Both files contain the same finished plate, without separate centerline objects.

The SVG preview still plots centerlines, not the finished plate. At the
requested spacing, the nominal line pitch exceeds the tool diameter by
0.0125 inch, leaving ridges between passes. Inspect the OBJ for the cut result.

For the original workflow, use `--geometry centerlines`: it exports
`centerlines.step`, with open spline edges and no solids. Adding `--with-plate`
in that mode also exports a separate `plate_reference.step`.

Generated outputs and the virtual environment are ignored by Git. Rerunning the
same output directory replaces its generated files. A failed validation/export
leaves existing artifacts intact. `report.json` lists the files produced in the
current run. Obsolete geometry files listed in the previous report are removed
after a successful run; unrelated files in the folder are preserved.

## Code layout and verification

- `src/millcraft/sampling.py`: image loading, coordinate grid, and depth samples.
- `src/millcraft/splines.py`: interpolation and analytical Z validation.
- `src/millcraft/step.py`: exact CAD curve translation and plate/STEP export.
- `src/millcraft/machining.py`: rectangular spline pockets, optional cutter envelopes, and OBJ export.
- `src/millcraft/preview.py`: portable SVG depth plot.
- `src/millcraft/cli.py`: workflow, CSV, and validation report.

```sh
.venv/bin/python -m pytest -q
```

Tests cover bilinear mapping, image orientation, transparency, aspect ratio,
spacing and input validation, an adversarial cubic that exceeds maximum depth
between sample knots, local overshoot within the global depth range, and a
safe natural cubic. The supplied-A integration test reimports the complete STEP
and verifies 41 B-spline edges, millimeter units, every sample and analytical
extremum, absence of solids, and the separate plate dimensions and volume in
centerlines mode. Pocket tests check rectangular outlines, straight ends,
identical floor depths across the width, exact spline-integral removal volumes,
shallow cuts, overlapping pockets, depth limits, and minimal face counts. Export
tests verify STEP round trips, closed oriented OBJ meshes, and failed-run
preservation. The optional cutter model retains capsule/rounded-end tests.
