#!/usr/bin/env bash
set -euo pipefail

# Run from the model directory, regardless of where this script is launched.
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

# Edit these settings, then run: bash run_2d.sh
PYTHON=python3
THREADS=2

# Domain and electrode (meters).
DOMAIN_R=2.0
DOMAIN_Z=4.0
COARSE_NR=16
COARSE_NZ=32
ROD_R0=0.0
ROD_Z0=0.0
ROD_R1=0.0
ROD_Z1=0.05
ROD_RADIUS=5e-3
VOLTAGE=-4e6

# Electric-field mesh refinement (spacing in meters, field in V/m).
MIN_DX=2e-4
MAX_DX=0.2
ELECTRODE_DX=1e-4
REFINE_E=3e6
DEREFINE_E=2e6
DEREFINE_LEVELS=2

# Permanent fine region beyond the physical rounded electrode tip.
# WIDTH is the full diameter: 0.002 means radial extent 0.001 m.
TIP_LENGTH=0.01
TIP_WIDTH=0.002
TIP_DX=1e-4

# Gas refinement outside the fixed region.
GAS_DX=2.5e-4
GAS_TEMPERATURE=301
GAS_AXIS_RADIUS=0        # 0 disables the full-length fine cylinder.
GAS_VARIATION=0.05
GAS_DYNAMICS=true

# Time integration and output.
DT=0.3e-9               # For a timestep comparison, try 0.15e-9.
N_STEPS=8000            # For a quick check, try 3 or 60.
STEPS_PER_OUTPUT=30
MEMORY_GB=64
POISSON_ATOL=10
VERBOSE=0               # Initialization progress is always displayed.
# Use a different output name for each comparison to avoid overwriting files.
OUTPUT_BASE=./output/output_2D/gas_4000kv_tip_focused_2d

export OMP_NUM_THREADS="$THREADS"

args=(
    -adaptive_grid
    -domain_size "$DOMAIN_R" "$DOMAIN_Z"
    -coarse_grid_size "$COARSE_NR" "$COARSE_NZ"
    -min_dx "$MIN_DX"
    -max_dx "$MAX_DX"
    -max_dx_electrode "$ELECTRODE_DX"
    -refine_E "$REFINE_E"
    -derefine_E "$DEREFINE_E"
    -derefine_nlevels "$DEREFINE_LEVELS"
    -tip_refine_length "$TIP_LENGTH"
    -tip_refine_width "$TIP_WIDTH"
    -tip_refine_dx "$TIP_DX"
    -gas_refine_dx "$GAS_DX"
    -gas_refine_temperature "$GAS_TEMPERATURE"
    -gas_refine_radius "$GAS_AXIS_RADIUS"
    -gas_refine_rtol "$GAS_VARIATION"
    -rod_r0 "$ROD_R0" "$ROD_Z0"
    -rod_r1 "$ROD_R1" "$ROD_Z1"
    -rod_radius "$ROD_RADIUS"
    "-phi_bc=$VOLTAGE"
    -poisson_atol "$POISSON_ATOL"
    -dt "$DT"
    -n_steps "$N_STEPS"
    -steps_per_output "$STEPS_PER_OUTPUT"
    -memory_limit "$MEMORY_GB"
    -verbose "$VERBOSE"
    -siloname "$OUTPUT_BASE"
)

case "$GAS_DYNAMICS" in
    true) args+=(-gas_dynamics) ;;
    false)
        if [[ "$GAS_DX" != 0 ]]; then
            echo 'Set GAS_DX=0, GAS_TEMPERATURE=0 and GAS_AXIS_RADIUS=0 when disabling gas dynamics.' >&2
            exit 2
        fi
        ;;
    *) echo 'GAS_DYNAMICS must be true or false.' >&2; exit 2 ;;
esac

# Additional command-line options override the settings above, for example:
# bash run_2d.sh -n_steps 3 -siloname /tmp/my_2d_check/run
exec "$PYTHON" -u 2d_model.py "${args[@]}" "$@"
