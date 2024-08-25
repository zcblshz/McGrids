# McGrids: Monte Carlo-Driven Adaptive Grids for Iso-Surface Extraction

This repository contains the official implementation of the ECCV paper **"McGrids: Monte Carlo-Driven Adaptive Grids for Iso-Surface Extraction"**. McGrids introduces a novel approach to enhance the efficiency of iso-surface extraction by constructing adaptive grids, as opposed to the traditional uniform grids used in prior work.

## Overview

McGrids leverages a Monte Carlo process to solve the problem of constructing adaptive grids as a probability sampling problem. The result is a significant reduction in the number of implicit field queries, which leads to substantial memory savings while producing high-quality meshes with detailed geometry. We validate McGrids through extensive experiments, including both analytical signed distance functions (SDFs) from surface meshes and learned implicit fields from real multiview images.

## Dependencies

Before proceeding with the installation, ensure that you have installed the following dependency:

- [Geogram](https://github.com/BrunoLevy/geogram): Please follow the instructions provided in the Geogram repository to install it. After installation, update the path to Geogram in `./differentiable_mcmt/CMakeLists.txt`.

## Installation

To set up the environment and install McGrids, follow these steps:

```bash
conda create -n mcgrids python=3.10.0
conda activate mcgrids
pip install .
```

To extract a mesh from an SDF using McGrids, you can run the following example:
```bash
cd examples
python extract_mesh_from_sdf.py --input_path ../assets/armadillo.obj --output_path ../assets/armadillo_mcgrids.obj
