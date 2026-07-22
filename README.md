# Geometric Digital Twinning Protocol for axisymmetric shells of revolution - Version 1, July 2026
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
3) Run protocolRun.m to generate a geometric digital twin of the sample WTST. A standard diagnostic for every stage is saved as an interactive MATLAB `.fig` file and a 300-dpi `.png` preview in the same input directory as the point-cloud and cached MAT files. Axes, tick labels, legends and colour bars use MATLAB's LaTeX interpreter.
4) Run algorithms in the ToleranceAssessment directory to assess buckling-relevant fabrication tolerances of the sample WTST to EN 1993-1-6 (2025).

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

# Reproducibility and cached stages
`protocolRun.m` defines a single random seed and stores the complete stage
configuration and point-cloud file signature in each cached MAT-file. A cache
is reused only when this metadata matches the current run. Increment
`implementationVersion` in `protocolRun.m` after code changes that should
invalidate existing results. Detailed transient plotting within individual
stages is disabled by default to avoid substantial run-time and memory
overhead; set `plotting = true` when those figures are required. Standard
summary diagnostics for all six stages are controlled separately by `generateDiagnostics` and
are regenerated only when missing or older than their corresponding MAT file.
Set `forceDiagnostics = true` to rebuild them unconditionally. They can also
be rebuilt without rerunning the protocol, for example:
```matlab
generateProtocolDiagnostics('ReleaseTower',8,53,true)
```

# Links
- Sample WTST dataset a.k.a. the 'release tower': https://doi.org/10.6084/m9.figshare.28903136.v1

# Publications
- Kathirkamanathan L, Sadowski AJ and Seidel M (2026) "Geometric digital twinning of full-scale metal shells of revolution - modelling protocol and buckling-relevant fabrication tolerance assessment" Under review.