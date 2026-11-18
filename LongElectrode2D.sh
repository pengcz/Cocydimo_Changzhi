#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

# Match the Python 3.14 extension; override only after rebuilding for another Python.
exec "${COCYDIMO_PYTHON:-/usr/bin/python}" 2d_model.py \
  -L_E_min 0.00005 \
  -L_E_max 0.025 \
  -domain_size 2 4 \
  -rod_r0 0 0 \
  -rod_r1 0 1 \
  -rod_radius 0.005 \
  -phi_bc=-1.2e6 \
  -dt 3e-10 \
  -n_steps 4000 \
  -min_dx 0.0002 \
  -max_dx 0.2 \
  -max_dx_electrode 0.0002 \
  -refine_E 3.0e6 \
  -verbose 0 \
  -steps_per_output 30 \
  -memory_limit 16 \
  -siloname ./output/LongElectrode_2D_ \
  -rng_seed 42 \
  -adaptive_grid \
  "$@"
