# Geometric Digital Twinning Protocol for axisymmetric shells of revolution - Version 1, September 2026
Geometric digital twinning protocol for full-scale metal shells of revolution and algorithms for assessing the buckling-relevant fabrication tolerances to EN 1993-1-6 (2025).

# Authors
Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski  
Department of Civil and Environmental Engineering, Imperial College London, UK 

Dr Marc Seidel  
Siemens Gamesa Renewable Energy, Germany

# Description
This repository presents a standardised protocol for the systematic processing of reality capture surveys of shells of revolution in point cloud format. The protocol performs enhanced registration, can segmentation, systematic and random outlier removal, surface reconstruction, surface inpainting and projection to arbitrary meshes. The protocol facilitates the post-construction quality assessment, damage and lifetime extension evaluations, and advanced structural analysis with real geometries of shells of revolution. Novel algorithms for automatically assessing EN 1993-1-6 (2025) buckling-relevant fabrication tolerances (dimple, out-of-roundness and unintended eccentricity) are provided. The protocol is illustrated on a real high-definition wind turbine support tower (WTST) laser scan and is readily generalisable to other metal shells of revolution. A link to a 5.7 GB augmented sandbox WTST point cloud dataset to be used in conjunction with this repository can be found below.

# Requirements
Both the geometric digital twinning protocol and tolerance algorithms were developed and tested using MATLAB R2023a. The following MATLAB Add-Ons are required to be installed:
- Computer Vision Toolbox
- Global Optimization Toolbox
- Image Processing Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox

# Usage
1) Clone repository
2) Download the WTST dataset from https://doi.org/10.6084/m9.figshare.28903136.v1 into the Input_ReleaseTower directory and unzip.
3) Review the stage-labelled choices in the initial `USER SETTINGS` section of `protocolRun.m`, then run it to generate a geometric digital twin of the sample WTST. A standard diagnostic for every stage is saved as an interactive MATLAB `.fig` file and a 300-dpi `.png` preview in the same input directory as the point-cloud and cached MAT files. Axes, tick labels, legends and colour bars use MATLAB's LaTeX interpreter.
4) Run algorithms in the ToleranceAssessment directory to assess buckling-relevant fabrication tolerances of the sample WTST to EN 1993-1-6 (2025).

# Height range
`protocolRun.m` defaults to `minZ = -Inf; maxZ = Inf`. Before processing,
these limits are replaced with the lowest and highest measured elevations
across all input `.bin` files. You can also leave just one end open, or supply
two finite limits such as `minZ = 8; maxZ = 53` to crop the data.

The first automatic lookup reads the scans in small chunks and saves
`S0_PointCloudHeight.mat` beside them. Later lookups reuse it while the input
file signature is unchanged. The resolved finite limits are included in stage
cache checks, so changing from a cropped range to full height rebuilds the
affected stages. This uses the measured extent, including any outlying points;
it does not infer missing data. Stage 2 derives its joint count from the
internal nominal boundaries in `S0_Data_<tower>.m` within this range. In that
file, `IS_FLANGE` must contain one logical (or 0/1) value per nominal segment,
in the same order as `Z_BOTS`. `NAME` is only a label; names and thicknesses
are not used to infer segment types. For the release tower, entries 8 and 9
are true and the other entries are false. Missing or invalid flags cause an
error before registration starts. The shared boundary of two adjacent flange
segments is a flange-to-flange joint;
a boundary between a flange and a shell segment is treated as a weld.
Flange joints are matched to detected peaks by elevation, with a check for
their adjacent welds. Missing or ambiguous flange matches stop segmentation.
`cloudStrakes` is derived automatically from nominal geometry and the height
range. Only segments with both internal bounding junctions inside that range
are eligible. Tower-end segments and segments cut by a height limit are excluded;
their missing boundaries are not invented. For the full-height release tower,
this selects nominal segments 2--14; end segments 1 and 15 are excluded.
Missing interior junctions stop validation. Internal joint bounds retain their nominal
indices when passed to reconstruction and diagnostics.

Stage 2 keeps the existing peak calculation and selects the strongest candidate
within `jointSearchRadius` (metres) of each nominal junction. Search windows
are clipped at adjacent nominal midpoints so a peak cannot serve two joints.
A window with no candidate returns `NaN`; no junction is
invented. The saved `jointMatches` table records search bounds, candidate counts,
selected elevations, offsets and peak strengths. Validation also checks the
selected elevations using `jointNominalTolerance` (metres), clipped at the same
midpoints. Before Stage 3, the existing strake assignment must also agree with
the nominal bounds. Missing bounding joints, including unhandled tower ends,
stop the run. These checks apply to cached results as well.

With Stage 2's `trimFailedEndJoints = true`, missing junctions at either end
are excluded and only strakes bounded by the remaining contiguous junctions
are reconstructed. Missing interior junctions still stop validation. Set the
flag to `false` to require every junction. Saved arrays retain their nominal
rows, with `NaN` bounds for excluded ends and a `retainedJoints` mask. The
Stage 2 diagnostic shades the region outside the retained joint span.

# Building the point-cloud converter
`Common/PTS2BIN.exe` must be a native executable for the operating system on
which MATLAB is running. From the `Common` directory, compile it as follows.

Windows with 64-bit MinGW-w64:
```text
g++ -O2 -std=c++17 -static -static-libgcc -static-libstdc++ PTS2BIN.cpp -o PTS2BIN.exe
```

Linux with GCC:
```text
g++ -O2 -std=c++17 PTS2BIN.cpp -o PTS2BIN.exe
```

The `.exe` name is retained on Linux because `convertPTS2BIN.m` invokes that
filename, although the generated file is a native Linux executable.

The converter accepts whitespace-separated X, Y, Z coordinates, optionally
preceded by a single unsigned integer giving the point count. Extra columns
(such as intensity, colour or normals) are ignored. Blank lines and an initial
UTF-8 byte-order mark are accepted. If a count is provided, it must match the
number of coordinate rows; malformed rows or count mismatches fail conversion
and remove the incomplete output.

# Reproducibility and cached stages
`protocolRun.m` defines a single random seed and stores the complete stage
configuration and point-cloud file signature in each cached MAT-file. A cache
is reused only when this metadata matches the current run. Increment
`implementationVersion` in `protocolRun.m` after code changes that should
invalidate all stage caches, or `reconstructionVersion` for changes limited to
Stages 3--6. Use `startFromStage` to resume with compatible earlier results;
Stages 3--5 run together. Detailed transient plotting within individual
stages is disabled by default to avoid substantial run-time and memory
overhead; set `plotting = true` when those figures are required. Standard
summary diagnostics for all six stages are controlled separately by `generateDiagnostics` and
are regenerated only when missing or older than their corresponding MAT file.
With `generateDiagnostics = true` (the default), Stage 5 automatically saves
its surface and point-cloud diagnostics before Stage 6 starts, including
`DTP_Diagnostic_Stage5_PointCloud_<tower>.fig` and `.png`,
showing registered measured points in 3D and unwrapped circumferential coverage,
coloured by scan filename. This deterministic display sample contains up to
50,000 points per scan before outlier removal, so gaps remain visible alongside
the reconstructed-surface diagnostic. It is a visual coverage check, not an
automatic connectivity test.
All six stage diagnostics include short captions explaining the plotted data
and what to check. Figures are also refreshed when the diagnostic code changes.
Stage 2 saves `S2_JointDetection_<tower>.mat` and plots detected versus nominal
junctions before validation, when `generateDiagnostics` is enabled. A rejected
detection therefore remains available in the Stage 2 FIG/PNG, clearly labelled
as unvalidated, without overwriting the completed segmentation cache. Once
segmentation succeeds, the plot is updated to include the fitted joint bounds.

Stage 4 initialises its local cone fit at the aligned tower axis and checks
convergence and positive radii. It uses a fixed neighbour-search radius of
`searchRadiusGridSpacingSF` times the estimated point spacing. Grid nodes with
no measured neighbours within that radius remain `NaN`; the radius does not
expand across missing regions. Smoothing preserves this missing-data mask,
and missing areas propagate through Stage 5 and mesh projection. The resulting
mesh can therefore be incomplete. Reconstruction caches from the earlier
expanding-radius implementation are invalidated automatically on the next run.
Set `forceDiagnostics = true` to rebuild them unconditionally. They can also
be rebuilt without rerunning the protocol, for example:
```matlab
generateProtocolDiagnostics('ReleaseTower',8,53,true)
```

# Surface outputs
`S5_Surface_<tower>.mat` contains the reconstructed outer surface by strake,
nominal geometry, thicknesses and retained strake indices. It is sufficient
for an analyst to construct their own mesh; `S6_Mesh_<tower>.mat` contains the
resampled grid. Lengths in these outputs are in millimetres and angles in
radians. Unsupported regions remain `NaN`. Neither file is a complete FE
model: element connectivity, shell reference-surface treatment, materials and
boundary conditions must be supplied separately.

# Links
- Sample WTST dataset a.k.a. the 'release tower': https://doi.org/10.6084/m9.figshare.28903136.v1

# Publications
- Kathirkamanathan L, Sadowski AJ and Seidel M (2026) "Geometric digital twinning of full-scale metal shells of revolution - modelling protocol and buckling-relevant fabrication tolerance assessment" Under review.

# AI-assisted development
OpenAI Codex was used throughout the development and refinement of this
codebase, including code generation, review, debugging, testing and
documentation. Responsibility for the scientific methodology, released code
and its scientific use remains with the authors.
