---
oneliner: MATLAB implementation of the Rectangular Surface Parameterization paper (Corman & Crane, TOG 2025)
tags: [surface-parameterization, quad-meshing, cross-field, frame-field, matlab, octave, geometry-processing, quantization, seamless-maps]
stack: [MATLAB, GNU Octave, C++, CMake]
generated: 2026-09-06
commit: 5b93eee
placeholder: true
---
Reference implementation for the SIGGRAPH/TOG paper "Rectangular Surface Parametrization" by Etienne Corman and Keenan Crane. Computes cross/frame fields (trivial, curvature-aligned or smooth), a seamless parametrization with configurable distortion/Chebyshev/alignment energies, and an optional integer quantization step (via a bundled C++ tool by Coudert-Osmont et al.) so the result can be turned into a quad mesh. Working research code with example meshes and an Octave compatibility layer plus validation tests added on top of the paper's original MATLAB source.
