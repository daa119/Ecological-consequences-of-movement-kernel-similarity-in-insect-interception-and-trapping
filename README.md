# Ecological consequences of movement kernel similarity in insect interception and trapping

**MATLAB simulation code accompanying the research manuscript by Danish A. Ahmed**

## Overview

This repository contains the MATLAB code used to generate Figures 1–8 of the manuscript *Ecological consequences of movement kernel similarity in insect interception and trapping*.

The study asks whether mathematically similar short- and long-distance movement kernels produce similar ecological outcomes. Thin-tailed and heavy-tailed step-length distributions are linked through a common survival probability, *p*, and compared using an optimal relationship obtained by minimising the squared Euclidean distance between their probability densities.

The resulting comparisons are examined in three individual-based insect-movement models:

- **One-dimensional trap cropping:** Pest redistribution, trap-crop occupancy, emigration and removal.
- **Two-dimensional pitfall trapping:** Capture of ground-dwelling arthropods by different trap arrangements.
- **Three-dimensional flight interception:** Contact with a finite trapping panel, capture and edge leakage.

The central conclusion is that **optimal mathematical similarity between movement kernels does not guarantee ecological equivalence**. Differences in displacement structure interact with attraction, boundaries, trap arrangement and target geometry.

## Repository contents

The repository contains seven self-contained MATLAB files corresponding to the eight manuscript figures.

| File | Manuscript figure | Purpose |
|---|---|---|
| `figure_1.m` | Figure 1 | Spatial distributions of insects at 0, 10 and 20 days under Gaussian and optimally related Cauchy movement in a trap-cropping system. |
| `figure_2.m` | Figure 2 | Time courses of cash-crop and trap-crop occupancy, cumulative removal and emigration for Gaussian movement and Cauchy movement at *p* = 0.1, 0.5, 0.721 and 0.9. |
| `figure_3.m` | Figure 3 | Heat maps of the difference in 20-day cumulative removal (long-distance minus short-distance) as trap-crop allocation and attraction strength vary. |
| `figure_4.m` | Figure 4 | Ground-arthropod spatial snapshots at 0, 30 and 60 days in a circular domain, comparing 3 × 3 grid and nine-trap transect arrangements. |
| `figure_5.m` | Figure 5 | Mean cumulative pitfall-trap capture over 30 days, with pointwise 95% confidence intervals from 10 simulation realisations. |
| `figure_6_and_7.m` | Figures 6 and 7 | Three-dimensional insect flight trajectories, representative panel interceptions and the corresponding two-dimensional projections. |
| `figure_8.m` | Figure 8 | Effects of square-panel side length on capture at two sticking probabilities, panel interception and edge leakage. |

## Modelling framework

### 1D: Trap cropping (Figures 1–3)

Movement within a 100 m field is represented by Gaussian and Cauchy displacements. Attraction increases spatially towards the trap crop; insects entering the trap crop switch to a reduced-scale, unbiased Gaussian process. Removal can occur on entry or during residence, while individuals outside the field can return.

The optimal one-dimensional kernel relationship is **p* = 0.721**.

### 2D: Pitfall trapping (Figures 4–5)

Individuals move within a confined circular domain containing nine pitfall traps arranged either in a grid or along a transect. Gaussian displacement components give Rayleigh step lengths, while long-distance movement uses folded-Cauchy step lengths.

Trap encounters are checked along each movement segment using 100 intermediate sample positions.

The optimal two-dimensional relationship is **p* = 0.658**.

### 3D: Flight interception (Figures 6–8)

Flying insects approach a finite square trapping panel. Short-distance step lengths follow a Chi distribution with three degrees of freedom; long-distance step lengths follow a folded-Cauchy distribution.

Directions before the first panel-plane crossing are sampled using a target-oriented von Mises–Fisher distribution. Segment–panel intersections, probabilistic retention after contact, repeated contacts and crossings outside the finite panel are represented explicitly.

The optimal three-dimensional relationship is **p* = 0.650**.

## Requirements and running the code

- **Software:** MATLAB. The source files for Figures 3, 6–7 and 8 explicitly indicate MATLAB R2016a compatibility.
- **Inputs:** No external datasets or separately supplied helper files are referenced; the required helper functions are included within the relevant files.
- **Outputs:** Running each file generates the corresponding figure or figures in MATLAB. Figures are displayed, not automatically saved to image files.
- **Randomness:** Fixed random seeds are specified for reproducibility (`rng(1)` in Figures 1–4 and 6–8; `rng(25)` in Figure 5). Results are stochastic and may depend on MATLAB version and execution environment.

Download or clone the repository, set the MATLAB Current Folder to the repository directory, and run the file corresponding to the desired manuscript figure.

The commands are:

```matlab
figure_1
figure_2
figure_3
figure_4
figure_5
figure_6_and_7
figure_8
```

Each file can be run independently.

**Figure 3's parameter sweep and the repeated-trajectory simulations may be computationally intensive.** Model parameters and random seeds are specified near the beginning of each file.

## Interpretation

The simulations illustrate three context-dependent outcomes:

1. **Trap cropping:** Long-distance movement generally increases removal under weak attraction, but short-distance movement can become more effective as attraction strengthens.
2. **Pitfall trapping:** Long-distance movement generally increases cumulative capture in the confined domain, although some non-optimal parameter relationships can produce capture similar to short-distance movement.
3. **Flight interception:** Long-distance movement reduces panel interception and capture while increasing the frequency of crossing the panel plane outside the finite target.

The codes are provided to support transparency, reproducibility and further investigation of how movement-kernel assumptions influence ecological predictions.

## Citation

If you use or adapt the code, please cite the accompanying manuscript and the archived software release when available:

Ahmed, D. A. (2026). *Ecological consequences of movement kernel similarity in insect interception and trapping*. Manuscript.

**Zenodo software DOI:** To be added after the first archived GitHub release.

## Author

**Danish A. Ahmed**

Department of Mathematics and Natural Sciences  
Gulf University for Science and Technology  
Hawally, Kuwait

Email: Ahmed.D@gust.edu.kw

## License

See the repository's `LICENSE` file for the applicable terms of use, once added.
