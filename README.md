# Descan Correction for 4D-STEM (NBD / SPED) Data

A Jupyter notebook workflow for correcting **descan error** in 4D-STEM
datasets acquired in NBD (nanobeam diffraction) or SPED (scanning precession
electron diffraction) mode.

## Background

In 4D-STEM, the electron beam is scanned across the sample while a full
diffraction pattern (DP) is recorded at every scan (probe) position. Ideally
the direct beam / diffraction disk should land on the **same detector
pixel** at every probe position, since diffraction geometry shouldn't depend
on where the probe sits on the sample. In practice, imperfect alignment of
the scan coils and descan coils causes the direct-beam position to drift
smoothly across the detector as a function of scan position — usually a
close-to-linear "tilt" across the field of view. This notebook measures and
removes that drift.

## Why sample-and-fit instead of correcting every pattern independently?

Tools such as [`pyxem`'s `Diffraction2D.center_direct_beam`](https://pyxem.readthedocs.io/en/stable/reference/generated/pyxem.signals.Diffraction2D.center_direct_beam.html)
estimate the beam position **independently, pattern by pattern** (via
`cross_correlate`, `blur`, `interpolate`, or `center_of_mass`) and shift each
pattern to the detector center. That works well when every pattern has a
clean, unambiguous direct beam — but for crystalline NBD/SPED data this
often breaks down:

- Strong Bragg disks can be as bright as, overlap with, or outshine the
  direct beam, so per-pattern estimators can lock onto the wrong disk.
- Vacuum, amorphous, or strongly diffracting regions may have a weak,
  saturated, or missing direct beam.
- Per-pattern noise turns a physically smooth descan drift into a
  spatially-incoherent, per-pixel "correction."

Descan error is a property of the scan/descan coil alignment, **not of the
sample**, so it is a smooth, close-to-linear function of scan position. This
notebook exploits that by:

1. Letting the user hand-pick a modest number of scan positions where the
   beam center is unambiguous (e.g. vacuum, or a low-index zone axis with a
   clear central disk).
2. Statistically rejecting any of those points that disagree with the rest.
3. Fitting a single smooth plane through the survivors and evaluating it at
   *every* scan position.

The result is a descan-correction map that only needs a representative
subset of patterns to have a clean beam, rather than every single one — more
robust for beam-sensitive or strongly diffracting crystalline samples, at
the cost of some manual point selection and the assumption that descan error
is well described by a plane (true for most misalignments, but not for
higher-order/non-linear distortions). See the [Comparison with
`pyxem`](#comparison-with-pyxem) section below for when the fully-automatic
`pyxem` route is the better choice.

## Workflow

The notebook proceeds in eight steps:

| Step | What happens |
|---|---|
| 0 | Load the 4D dataset (`.hspy`) and build a real-space overview image (mean over the DP dimensions). |
| 1 | Interactively sample the direct-beam center at a sparse set of scan positions, using **one** of three methods (see below). |
| 2 | Visualize the sparse `descan_x` / `descan_y` maps (unsampled positions are `NaN`). |
| 3 | Robustly reject outlier points (mis-clicks, bad patterns, failed detections) via `remove_plane_outliers`. |
| 4 | Fit a plane `z = c0 + c1*x + c2*y` through the cleaned points with `createSurfaceFits`, and evaluate it at every scan position to get a full-field descan map (`descan_bg_x`, `descan_bg_y`). |
| 5 | Visualize the fitted, full-field descan planes. |
| 6 | Shift every diffraction pattern so its beam lands at the global mean beam position. |
| 7 | Save the corrected 4D dataset as `.hspy`, re-attaching the original calibration/metadata. |
| 8 | Verify the correction: re-inspect DPs interactively, and compare the PACBED pattern before/after (it should sharpen and re-center). |

### Step 1 — three interactive center-finding methods

Only run **one** of these per correction pass — each (re)initializes
`descan_x` / `descan_y`. All three share the same interaction: click a scan
position on the overview image (left panel) to load its DP into the right
panel, then interact with the DP to record a beam center.

- **Method A — Center of Mass (COM) in a hand-drawn ROI.** Draw a rectangle
  around the direct beam; the intensity-weighted center of mass is computed
  only inside that rectangle, excluding bright Bragg spots or contamination
  elsewhere in the pattern. Most robust when the direct beam is a diffuse
  blob without sharp edges, since COM doesn't need a well-defined boundary.
- **Method B — Direct manual click.** Click directly on the beam center.
  Fastest and most transparent (what you click is what you get) — useful as
  a sanity check against the other two methods — but limited to pixel-level
  precision, not sub-pixel.
- **Method C — Automatic circle detection (Canny edges + Hough transform).**
  Draw a rough ROI around the direct beam; edges are detected (Canny) and a
  circle is fit via the Hough circle transform. The most automated and most
  precise (sub-pixel, symmetric) of the three when the beam has a clean
  circular edge, but sensitive to the expected radius range (`min_r`,
  `max_r`) and can fail silently ("No circle detected") if the beam is
  defocused, elliptical, or overlaps other disks.

You only need to sample a modest number of positions (tens, not thousands),
spread across the whole scan area, at positions where the direct beam is
unambiguous.

### Step 3 — outlier rejection (`remove_plane_outliers`)

1. Fits a first-pass plane to all sampled (non-NaN) points.
2. Computes each point's residual from that plane.
3. Converts residuals to a robust z-score using the **median absolute
   deviation (MAD)** instead of the standard deviation, so a few bad points
   don't skew the threshold.
4. Discards (sets to `NaN`) any point whose robust z-score exceeds `thresh`
   (default `4.0`).

Applied independently to `descan_x` and `descan_y`.

### Step 4 — plane fit (`createSurfaceFits`)

Fits `z = c0 + c1*x + c2*y` by least squares through the outlier-cleaned
points (not the raw sampled points), and returns a callable that evaluates
the plane at any `(x, y)`, together with the fit coefficients and
goodness-of-fit metrics (R², adjusted R², RMSE). This plane is evaluated
over the full scan grid to produce `descan_bg_x` / `descan_bg_y`.

### Step 6 — applying the correction

For every scan position, the diffraction pattern is shifted (`scipy.ndimage.shift`,
linear interpolation, `nearest` boundary mode) by the difference between the
global mean beam position (mean of the raw sampled COMs) and the *fitted*
local beam position at that pixel — not the raw sampled value — so every
pattern ends up centered at the same reference point.

## Comparison with `pyxem`

| | This notebook | `pyxem.center_direct_beam` (default use) |
|---|---|---|
| Sampling | Sparse, user-selected scan positions | Every scan position, automatically |
| Center estimate per point | Chosen interactively (COM / manual / Hough) with visual confirmation | One automatic estimator applied uniformly |
| Outlier handling | Explicit robust (MAD-based) rejection before fitting | None built in |
| Model of descan error | Explicit smooth-plane assumption, enforced by least-squares fit | None by default — no cross-pattern smoothness constraint |
| Robustness on crystalline/strong-diffraction data | High — only needs a subset of patterns to have an unambiguous beam | Lower — every pattern needs the estimator to find the beam correctly |
| Effort | Manual point-picking required | Fully automatic, one function call |
| Speed | Slower (interactive + explicit Python loop) | Fast, vectorized/lazy-friendly |

`pyxem`'s `center_direct_beam` is preferable when the direct beam is
reliably the brightest, most well-defined feature in *every* pattern (e.g.
amorphous samples, diffuse/incoherent imaging), or when the true
beam-position drift is not well described by a plane. A hybrid is also
possible: `pyxem`'s `BeamShift.get_linear_plane` / `make_linear_plane`
methods fit a plane through *all* per-pattern estimates (optionally with
outlier masking) — the automated, dense-sampling analog of Steps 3–4 here,
letting you combine `pyxem`'s fast per-pattern estimation with the same
smoothness constraint used in this notebook. See the notebook's Appendix for
details.

## Requirements

- Python 3.11+
- [`hyperspy`](https://hyperspy.org/)
- `numpy`, `scipy`, `h5py`
- [`scikit-image`](https://scikit-image.org/) (`skimage.exposure`, `skimage.feature.canny`, `skimage.filters.gaussian`, `skimage.transform.hough_circle*`)
- `matplotlib` (with the `%matplotlib widget` / `ipympl` interactive backend)
- `tqdm`

Install with:

```bash
pip install hyperspy numpy scipy h5py scikit-image matplotlib ipympl tqdm
```

## Data format

Input is a 4D-STEM dataset stored as a HyperSpy `.hspy` file with shape
`(scan_x, scan_y, detector_x, detector_y)`.

## Usage

1. Open `Descan_Correction_NBD_SPED.ipynb` in Jupyter (JupyterLab or
   Notebook, with the `ipympl`/`widget` backend enabled for interactivity).
2. Update the `path` and `file_name` variables in the **Load the 4D-STEM
   Dataset** cell to point at your `.hspy` file.
3. Run the **Data Visualization** cells to inspect the scan overview and
   individual diffraction patterns.
4. Run **one** of the three Step 1 method cells and click through a spread
   of scan positions to sample the direct-beam center.
5. Run the Step 2–5 cells to visualize the sparse map, clean outliers, and
   fit/inspect the full-field descan plane.
6. Run Step 6 to apply the correction across the whole dataset.
7. Run Step 7 to save the corrected dataset (saved alongside the original as
   `<name>_dscan_corrected.hspy`).
8. Run Step 8 to verify the correction via the interactive DP viewer and the
   before/after PACBED comparison.

## Output

A descan-corrected 4D-STEM dataset saved as `<original_name>_dscan_corrected.hspy`,
with the original axes/calibration and metadata preserved.

## License

Add a license of your choice (e.g. MIT, BSD-3-Clause) here.
