cocydimo: Conducting Cylinder Discharge Model
=======

This code implements a reduced discharge model for streamer discharges, in which the channels are represented by a collection of conducting cylindrical segments. The conductivity and electric potential together evolve on a numerical mesh, as described in this [paper](https://doi.org/10.1016/j.cpc.2025.109733).

Getting the code and compiling
==

To compile the code, the following packages are required:

* Recent versions of `gcc` and `gfortran`. Other compilers can probably be used as well by modifying the build process.
* `f2py` and `meson` and python3 development files (e.g., `python3-dev' on Debian/Ubuntu) to create a Python module

The code (including the [afivo library](https://github.com/MD-CWI/afivo)) can be downloaded with:

    git clone --recurse-submodules https://github.com/jannisteunissen/cocydimo.git

Compiling for the first time can take some time:

    cd cocydimo
    make

Running
==

The model can be run using

    ./2d_model.py # 2D axisymmetric, uniform grid
    ./2d_model.py -adaptive_grid # 2D axisymmetric, adaptive mesh refinement
    ./3d_model # 3D with adaptive mesh refinement and branching

To see a full list of options, use:

    ./2d_model.py -h
    ./3d_model.py -h

For 2D adaptive refinement, use `-min_dx`, `-max_dx`, and
`-max_dx_electrode` to control grid spacing, and `-refine_E` and
`-derefine_E` to control the electric-field thresholds. For example:

    python3 2d_model.py -adaptive_grid -min_dx 1e-4 -max_dx 2e-3 -max_dx_electrode 8e-4

The `-grid_size` option applies only to the default uniform grid. Adaptive
mode updates refinement each time step and refreshes the grid spacing used
for streamer sampling.

For gas heating studies in 2D, add:

    -gas_dynamics -gas_refine_dx 2.5e-4 -gas_refine_radius 0.02 -gas_refine_rtol 0.05

`-gas_refine_dx` enables gas refinement and sets its maximum cell spacing.
`-gas_refine_radius` keeps a cylinder around the axis fine along the entire
axial domain, with one box of radial buffer. This resolves the heating source
before temperature gradients develop and retains resolution after the tip
passes. Outside it, relative box variations in temperature, pressure, or
density above `-gas_refine_rtol` trigger refinement, including ghost cells
to detect features at box boundaries. Variation above half this threshold
prevents coarsening. Gas refinement overrides distance-from-tip coarsening.
Set the radius to zero for variation-only refinement, or the tolerance to
zero for cylinder-only refinement. Gas refinement is disabled by default.

To keep the whole warm region fine, including regions with almost uniform
temperature, add `-gas_refine_temperature 301`. Boxes containing gas above
301 K request `-gas_refine_dx` resolution, regardless of distance from the
tip. Ghost cells provide a one-cell buffer across box boundaries. This rule
retains existing refinement until temperature falls to 300.5 K or below;
other active criteria may still protect it. Set this option to zero to disable
the temperature threshold. Refinement changes by at most one level per time
step, so a newly heated coarse region takes multiple steps to reach its target.
Keep a modest permanent axis region to resolve heating before it starts.

The 2D driver always displays initialization refinement/Poisson progress and
flushes it immediately, including when redirected to a log. After initialization,
`-verbose` controls solver logging. Ctrl+C terminates native solver calls directly;
an interrupted output file may be incomplete. The 2D default `-poisson_atol` is
10 (previously 1): the 4 MV, approximately 0.061 mm configuration stalled at a
residual of 1.5 with the old absolute tolerance. Relative tolerance is unchanged.
This is a stopping-tolerance adjustment, not a temperature convergence check.

To focus resolution just beyond an axial electrode, use
`-tip_refine_length 0.01 -tip_refine_width 0.002 -tip_refine_dx 1e-4`.
This permanently resolves boxes intersecting the cylinder from the physical
rounded tip to 1 cm beyond it, with full diameter 2 mm (radial extent 1 mm).
The physical tip is at max(rod_r0_z, rod_r1_z) + rod_radius. Both electrode
endpoints must be on the axis. The fixed-region spacing is independent of
electric-field and gas spacings. Set `-gas_refine_radius 0` to avoid refining
the full axial domain. Temperature and variation criteria still operate outside
the fixed region if enabled. Block boundaries and mesh balancing can extend
the fine grid beyond the requested cylinder.

These controls require rebuilding the native solver. Use a fresh build
directory if an existing Meson directory was copied from another checkout:

    make lib2d BUILD_DIR_2D=build_gas_2d

For the 2 m by 4 m domain with coarse grid 16 by 32, a requested gas spacing
of 2.5e-4 m gives 0.244 mm cells. Halve the requested spacing to check
convergence of temperature peaks and profiles; this can substantially increase
runtime and memory. The permanent cylinder includes unheated gas ahead of
the streamer. Increase its radius if the hot region extends beyond it.

Example
==

![alt text](https://github.com/jannisteunissen/cocydimo/blob/main/media/ICPIG_2025_example.png?raw=true)

### Open lower boundary in 3D

Pass `-open_lower_boundary` to apply a local Robin condition directly at
the original lower face outside the grounded rod. The upper face uses the
applied voltage, the sides retain zero normal field, and the grid dimensions
are unchanged. Omitting the flag grounds the lower plane.

The approximation uses dφ/dz = k (φ − Vtop) at z = 0, with
k = π / max(Lx, Ly). It assumes the exterior potential tends towards the
applied top voltage and uses the decay rate of the longest nonconstant
lateral mode compatible with the insulating sides. A general open half-space
has different decay rates for different lateral modes. This local condition
therefore approximates exterior air and can affect quantitative tip fields.

```bash
./LongElectrode_3D.sh -open_lower_boundary -siloname ./output/LongElectrode_open_3D
```
