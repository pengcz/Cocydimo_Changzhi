#!/usr/bin/env python3
"""Cocydimo: Conducting Cylinder Discharge Model"""
#Reference Command line ./3d_model.py     -domain_size 4.0 4.0 4.0     -coarse_grid_size 32 32 32     -rod_r0 2.0 2.0 0.0     -rod_r1 2.0 2.0 0.05     -rod_radius 5e-3     -phi_bc=-200e4     -phi_factor_vs_time phi_factor.txt     -n_initial_streamers 1     -dt 0.3e-9     -n_steps 8000     -min_dx 2e-4     -max_dx 0.2     -max_dx_electrode 2e-3     -derefine_nlevels 2     -steps_per_output 30     -memory_limit 64     -refine_E 3.0e6     -verbose 0     -print_performance     -siloname /scratch-shared/cpeng/output_reduced_260924/output/2MV/gas_2000kv_dynamics_3d     -gas_dynamics

import copy
import argparse
import numpy as np
import json
import os
import signal
import sys
from time import perf_counter
from numpy.linalg import norm
import model_lib as mlib
# Python's deferred SIGINT handler cannot run during a long native solver call.
# Use the OS handler so Ctrl+C also stops initialization immediately.
signal.signal(signal.SIGINT, signal.SIG_DFL)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault('GFORTRAN_UNBUFFERED_PRECONNECTED', 'y')
from poisson_2d import m_solver as p2d

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Cocydimo: Conducting Cylinder Discharge Model')
parser.add_argument('-r_scale', type=float, default=1.2,
                    help='Scale factor compared to electrodynamic radius')
parser.add_argument('-domain_size', type=float, nargs=2,
                    default=[30e-3, 30e-3], help='Domain size (m)')
parser.add_argument('-coarse_grid_size', type=int, nargs=2,
                    default=[16, 16],
                    help='Size of coarse grid (Nr, Nz)')
parser.add_argument('-grid_size', type=int, nargs=2, default=[256, 256],
                    help='Size of computational grid (Nr, Nz)')
parser.add_argument('-adaptive_grid', action='store_true',
                    help='Use adaptive refinement instead of the uniform grid')
parser.add_argument('-tip_refine_length', type=float, default=0.0,
                    help='Length beyond the physical electrode tip to keep fine (m); 0 disables')
parser.add_argument('-tip_refine_width', type=float, default=0.002,
                    help='Full diameter of the fixed fine region around the axis (m)')
parser.add_argument('-tip_refine_dx', type=float, default=1e-4,
                    help='Maximum cell spacing in the fixed tip region (m)')
parser.add_argument('-gas_refine_dx', type=float, default=0.0,
                    help='Maximum spacing in gas refinement regions (m); 0 disables')
parser.add_argument('-gas_refine_radius', type=float, default=0.0,
                    help='Permanently resolve this radius around the axis (m), plus a box buffer')
parser.add_argument('-gas_refine_rtol', type=float, default=0.05,
                    help='Relative box variation in temperature, pressure or density triggering gas refinement')
parser.add_argument('-gas_refine_temperature', type=float, default=0.0,
                    help='Refine gas hotter than this temperature (K); 0 disables; retain to 0.5 K below')
parser.add_argument('-refine_E', type=float, default=3e6,
                    help='Refine above this electric field (V/m)')
parser.add_argument('-derefine_E', type=float, default=2e6,
                    help='Derefine below this electric field (V/m)')
parser.add_argument('-derefine_nlevels', type=int, default=1,
                    help='Derefine at most this many levels')
parser.add_argument('-min_dx', type=float, default=1e-4,
                    help='Minimum allowed grid spacing (m)')
parser.add_argument('-max_dx', type=float, default=2e-3,
                    help='Maximum allowed grid spacing (m)')
parser.add_argument('-max_dx_electrode', type=float, default=8e-4,
                    help='Maximum grid spacing around the electrode (m)')
parser.add_argument('-refine_max_dist_head', type=float, default=5e-3,
                    help='Maximum distance from an active head for refinement (m)')
parser.add_argument('-box_size', type=int, default=8,
                    help='Size of boxes in afivo (#cells)')
parser.add_argument('-rod_r0', type=float, nargs=2,
                    default=[0.0e-3, 0.0e-3],
                    help='First point of rod electrode (m)')
parser.add_argument('-rod_r1', type=float, nargs=2,
                    default=[0.0, 5.0e-3],
                    help='Second point of rod electrode (m)')
parser.add_argument('-rod_radius', type=float, default=0.75e-3,
                    help='Radius of rod electrode (m)')
parser.add_argument('-n_steps', type=int, default=10,
                    help='How many steps to simulate')
parser.add_argument('-dt', type=float, default=2.5e-10,
                    help='Time step (s)')
parser.add_argument('-dt_factor', type=float, default=1.0,
                    help='Increase dt by this factor when streamers are gone')
parser.add_argument('-dz_data', type=float, default=30e-3/256,
                    help='Grid spacing used to obtain L_E from dataset (m)')
parser.add_argument('-phi_bc', type=float, default=-4e4,
                    help='Applied potential (V)')
parser.add_argument('-phi_factor_vs_time', type=str,
                    help='File with factor for applied voltage vs time')
parser.add_argument('-use_circuit', action='store_true',
                    help='Use R-C circuit for voltage source')
parser.add_argument('-capacitance', type=float, default=1e-9,
                    help='Capacitance C of the voltage source (farad)')
parser.add_argument('-resistance', type=float, default=300.0,
                    help='Internal resistance R of the voltage source (Ohm)')
parser.add_argument('-alpha', type=float, default=0.5,
                    help='Exponential smoothing coefficient')
parser.add_argument('-channel_update_delay', type=float, default=1e-9,
                    help='Delay for first updating channel conductivity (s).'
                    'Should be larger than the dielectric relaxation time.')
parser.add_argument('-channel_no_ionization', action='store_true',
                    help='Do not increase channel conductivity, which can be '
                    'problematic near domain boundaries)')
parser.add_argument('-channel_max_sigma', type=float, default=5.0,
                    help='Limit growth of volume conductivity to this value '
                    'to prevent issues at domain boundaries [A/(m V)]')
parser.add_argument('-channel_min_sigma', type=float, default=1e-9,
                    help='Keep electron condtivity above this value to allow '
                    'it to grow again [A/(m V)]')
parser.add_argument('-mu_electron', type=float, default=0.04,
                    help='Effective electron mobility at 1 bar, 300 K. '
                    'mu_e/mu_i is used to update ion conductivity (m2/(V s))')
parser.add_argument('-mu_ion', type=float, default=2e-4,
                    help='Effective ion mobility at 1 bar, 300 K')
parser.add_argument('-k_ion_recombination', type=float, default=1e-13,
                    help='Ion-ion recombination rate constant (m^3/s)')
parser.add_argument('-L_E_max', type=float, default=5e-3,
                    help='Maximum of L_E (m) at 1 bar, 300 K (is scaled)')
parser.add_argument('-L_E_min', type=float, default=1e-4,
                    help='Minimum of L_E (m) at 1 bar, 300 K (is scaled)')
parser.add_argument('-c0_L_E_dx', type=float, default=0.75,
                    help='Correction factor for L_E w.r.t. data grid spacing')
parser.add_argument('-c1_L_E_dx', type=float, default=0.0,
                    help='Correction factor for L_E when dx < dx_data')
parser.add_argument('-transport_data_file', type=str,
                    default='data/TD_N2_0.8_Phelps_O2_0.2_Phelps.txt',
                    help='Transport data file')
parser.add_argument('-k_eff_num_points', type=int, default=200,
                    help='Number of points to use internally for k_eff_table')
parser.add_argument('-poisson_rtol', type=float, default=1e-5,
                    help='Relative tolerance for Poisson solver')
parser.add_argument('-poisson_atol', type=float, default=10.0,
                    help='Absolute tolerance for Poisson solver')
parser.add_argument('-siloname', type=str, default='output/simulation_2d',
                    help='Base filename for output Silo files')
parser.add_argument('-write_eps', action='store_true',
                    help='Write epsilon (of Poisson eq.) to output')
parser.add_argument('-write_time', action='store_true',
                    help='Write time that channel was added to output')
parser.add_argument('-write_rhs', action='store_true',
                    help='Write r.h.s. of Poisson eq. to output')
parser.add_argument('-rng_seed', type=int,
                    help='Seed for the random number generator')
parser.add_argument('-memory_limit', type=float, default=0.1,
                    help='Memory limit (GB)')
parser.add_argument('-steps_per_output', type=int, default=1,
                    help='Write output every N steps')
parser.add_argument('-gas_dynamics', action='store_true',
                    help='Simulate gas dynamics')
parser.add_argument('-pressure', type=float, default=1.0,
                    help='Gas pressure (bar)')
parser.add_argument('-temperature', type=float, default=300.0,
                    help='Gas temperature (Kelvin)')
parser.add_argument('-mean_molecular_weight', type=float, default=28.97,
                    help='Mean molecular weight of gas molecules (Dalton)')
parser.add_argument('-gas_gamma', type=float, default=1.4,
                    help='Gas adiabatic index')
parser.add_argument('-gas_fast_heat_factor', type=float, default=1.0,
                    help='Fraction of Joule heating that is immediately'
                    'converted to gas heating')
parser.add_argument('-gas_slow_heat_factor', type=float, default=0.0,
                    help='Fraction of Joule heating that is slowly converted'
                    'to gas heating')
parser.add_argument('-gas_slow_heat_timescale', type=float, default=20.0e-6,
                    help='Time scale for slow heating (s)')
parser.add_argument('-verbose', type=int, default=0,
                    help='How verbose the code is (> 0 shows more info)')

args = parser.parse_args()

k_b = 1.3806503e-23             # Boltzmann constant

if not np.isfinite(args.tip_refine_length) or args.tip_refine_length < 0:
    parser.error('-tip_refine_length must be finite and nonnegative')
if args.tip_refine_length > 0:
    if not args.adaptive_grid:
        parser.error('Tip region refinement requires -adaptive_grid')
    if not all(np.isfinite(x) and x > 0 for x in
               [args.tip_refine_width, args.tip_refine_dx]):
        parser.error('Tip region width and spacing must be finite and positive')
    if args.rod_r0[0] != 0 or args.rod_r1[0] != 0:
        parser.error('Tip region refinement requires an electrode on the axis')
    tip_z = max(args.rod_r0[1], args.rod_r1[1]) + args.rod_radius
    if tip_z < 0 or tip_z + args.tip_refine_length > args.domain_size[1]:
        parser.error('Tip refinement region must lie inside the axial domain')
    if args.tip_refine_width / 2 > args.domain_size[0]:
        parser.error('Tip refinement radius exceeds the radial domain')

gas_options = [args.gas_refine_dx, args.gas_refine_radius, args.gas_refine_rtol,
               args.gas_refine_temperature]
if not all(np.isfinite(x) and x >= 0 for x in gas_options):
    parser.error('Gas refinement parameters must be finite and nonnegative')
if args.gas_refine_dx > 0 and not (args.adaptive_grid and args.gas_dynamics):
    parser.error('-gas_refine_dx requires -adaptive_grid and -gas_dynamics')
if args.gas_refine_radius > 0 and args.gas_refine_dx == 0:
    parser.error('-gas_refine_radius requires a positive -gas_refine_dx')
if args.gas_refine_temperature > 0 and args.gas_refine_dx == 0:
    parser.error('-gas_refine_temperature requires a positive -gas_refine_dx')

if args.adaptive_grid:
    spacings = [args.min_dx, args.max_dx, args.max_dx_electrode]
    if not all(np.isfinite(dx) and dx > 0 for dx in spacings):
        parser.error('Refinement grid spacings must be finite and positive')
    # Electrode and gas criteria can request finer cells than the field criterion.
    if args.min_dx > args.max_dx:
        parser.error('-min_dx must not exceed -max_dx')
    if not (np.isfinite(args.refine_E) and np.isfinite(args.derefine_E)
            and 0 <= args.derefine_E <= args.refine_E):
        parser.error('Require 0 <= -derefine_E <= -refine_E (finite values)')
    if args.derefine_nlevels < 0:
        parser.error('-derefine_nlevels must be nonnegative')
    if not np.isfinite(args.refine_max_dist_head) or args.refine_max_dist_head <= 0:
        parser.error('-refine_max_dist_head must be finite and positive')

# Make sure output folder exists
os.makedirs(os.path.dirname(args.siloname), exist_ok=True)

# Save settings
fname = f'{args.siloname}.json'
with open(fname, 'w') as f:
    json.dump(args.__dict__, f, indent=2)
    print(f'Wrote settings to {fname}')

model = mlib.AirStreamerModel(c0=args.c0_L_E_dx, c1=args.c1_L_E_dx,
                              dz0=args.dz_data)


def set_voltage(time, V0, voltage_table):
    voltage = V0
    if voltage_table is not None:
        voltage *= np.interp(time, voltage_table[0], voltage_table[1])
    p2d.set_voltage(voltage)


if args.phi_factor_vs_time is not None:
    voltage_table = np.loadtxt(args.phi_factor_vs_time).T
else:
    voltage_table = None


np.random.seed(args.rng_seed)

def configure_solver_logging(verbosity):
    p2d.store_parameters(args.channel_min_sigma,
                         args.channel_max_sigma,
                         args.mu_electron,
                         args.mu_ion,
                         args.k_ion_recombination,
                         args.resistance,
                         args.capacitance,
                        verbosity)


# Show initialization progress even at the default simulation verbosity.
configure_solver_logging(max(1, args.verbose))
initialization_start = perf_counter()
print('Initializing domain and electric field; Ctrl+C stops the solver.')
print(f'Poisson tolerances: rtol={args.poisson_rtol:g}, atol={args.poisson_atol:g}')

p2d.set_rod_electrode(args.rod_r0, args.rod_r1, args.rod_radius)
p2d.initialize_domain(args.domain_size, args.coarse_grid_size,
                      args.box_size, args.phi_bc, args.memory_limit,
                      args.write_eps, args.write_time, args.write_rhs,
                      args.gas_dynamics)
if args.adaptive_grid:
    print('Building adaptive grid (refinement and Poisson iterations follow).')
    if args.tip_refine_length > 0:
        lower = [0.0, tip_z]
        upper = [args.tip_refine_width / 2, tip_z + args.tip_refine_length]
        print(f'Fixed fine region: r=0..{upper[0]:g} m, '
              f'z={lower[1]:g}..{upper[1]:g} m; dx <= {args.tip_refine_dx:g} m')
        p2d.set_focus_refinement(lower, upper, args.tip_refine_dx)
    if args.gas_refine_dx > 0:
        p2d.set_gas_refinement(args.gas_refine_dx, args.gas_refine_radius,
                              args.gas_refine_rtol)
        p2d.set_gas_temperature_refinement(args.gas_refine_temperature)
    p2d.set_refinement(args.refine_E, args.derefine_E,
                       args.min_dx, args.max_dx, args.max_dx_electrode,
                       args.derefine_nlevels, args.refine_max_dist_head,
                       args.poisson_rtol, args.poisson_atol)
else:
    print('Building uniform grid.')
    p2d.use_uniform_grid(args.grid_size)

print(f'Grid initialization finished in {perf_counter() - initialization_start:.1f} s.')
configure_solver_logging(args.verbose)

# Initial gas density
N0 = 1e5 * args.pressure / (args.temperature * k_b)

# Scaled values for L_E_max and min
L_E_max = args.L_E_max * (1.0/args.pressure) * (args.temperature/300)
L_E_min = args.L_E_min * (1.0/args.pressure) * (args.temperature/300)

if args.gas_dynamics:
    p2d.set_gas(args.pressure, args.temperature, args.mean_molecular_weight,
                args.gas_gamma, args.gas_fast_heat_factor,
                args.gas_slow_heat_factor, args.gas_slow_heat_timescale)

dz = p2d.get_finest_grid_spacing()
print(f'Minimum grid spacing: {dz:.2e}')

# Compute initial solution
set_voltage(0.0, args.phi_bc, voltage_table)
p2d.solve(0.0, args.poisson_rtol, args.poisson_atol)
p2d.write_solution(f'{args.siloname}_{0:04d}', 0, 0.)

# Set table with effective ionization rate
table_fld, table_k_eff = mlib.effective_ionization_rate(
    args.transport_data_file, args.temperature, args.pressure)

# Ensure table has uniform spacing
x = np.linspace(table_fld[0], table_fld[-1], args.k_eff_num_points)
y = np.interp(x, table_fld, table_k_eff)

if args.channel_no_ionization:
    y = np.minimum(y, 0.0)

p2d.store_k_eff(x[0], x[-1], y)

# Get L_E to estimate initial streamer radius
Emax, r_Emax = p2d.get_max_field_location()
z, E, success = p2d.get_var_along_line('E_norm', [0.0, r_Emax[1]], [0., 1.0],
                                       L_E_max, 2*L_E_max/dz)
if not success:
    raise RuntimeError('Interpolation error')

# Use the sampling length when the field stays above threshold along the line.
L_E = model.get_L_E(z, E, N0, dz, prev=L_E_max)

# Start with a smaller radius to approximate initial phase
radius0 = 0.5 * args.r_scale * model.get_radius(L_E, N0)

# Single streamer in the z-direction
streamers = [mlib.Streamer([0.0, r_Emax[1] - radius0],
                           [0., 1.0], radius0, 0.0)]
# Preserve the initial estimate for profiles that stay above threshold.
streamers[0].L_E = L_E

time = 0.0
V_cap = args.phi_bc
V_gap = args.phi_bc
step = 0

# Write current to a file
f_current = open(f'{args.siloname}_current.txt', 'w')
f_current.write('# time(s) J_tot J_displ V_gap V_cap\n')

J_tot, J_displ, G_eff = p2d.compute_current(time)
f_current.write(f'{time:.6e} {J_tot:.6e} {J_displ:.6e}\n')
f_current.flush()

for step in range(1, args.n_steps+1):
    print(f'{step:4d} t = {time*1e9:.1f} ns')

    # Potentially increase dt when there are no more streamers
    if len(streamers) == 0:
        dt = args.dt * args.dt_factor
    else:
        dt = args.dt

    streamers_prev = copy.deepcopy(streamers)

    for s in streamers:

        # Get samples of |E| ahead of the streamer to determine L_E
        r_tip = s.r + 0.5 * s.R * s.v/norm(s.v)

        # Check whether we can interpolate at r_tip
        E_hat, success = p2d.get_field_vector_at(r_tip)

        if not success:
            print('Could not sample E_hat, removing streamer')
            s.keep = False
            continue

        if args.gas_dynamics:
            N0, success = p2d.get_gas_number_density_at(r_tip)

        z, E, success = p2d.get_var_along_line('E_norm', r_tip, s.v,
                                               L_E_max, 2*L_E_max/dz)

        if success:
            L_E_new = model.get_L_E(z, E, N0, dz, prev=s.L_E)
        else:
            # Use previous value
            L_E_new = s.L_E

        if L_E_new < L_E_min:
            print(f'L_E too small {L_E_new:.2e}, removing streamer')
            s.keep = False
            continue

        if step == 1:
            L_E = L_E_new
        else:
            L_E = args.alpha * L_E_new + (1 - args.alpha) * L_E

        # Propagation in +z direction
        E_hat = np.array([0.0, 1.0])

        s.L_E = L_E
        s.sigma = model.get_sigma(L_E, N0)
        s.v = model.get_velocity(L_E, N0) * E_hat

        dR = min(args.r_scale * model.get_radius(L_E, N0) - s.R,
                 norm(s.v) * dt)
        s.R = s.R + dR
        s.r = s.r + s.v * (dt - 0.99 * dR/norm(s.v))

    if args.adaptive_grid:
        p2d.adjust_refinement()
        dz = p2d.get_finest_grid_spacing()

    if args.gas_dynamics:
        p2d.update_gas(dt)

    mlib.update_sigma(2, p2d.update_sigma, streamers, streamers_prev,
                      time, dt, args.channel_update_delay, step == 1)
    time += dt
    set_voltage(time, args.phi_bc, voltage_table)
    p2d.solve(dt, args.poisson_rtol, args.poisson_atol)

    J_tot, J_displ, G_eff = p2d.compute_current(time)

    if args.use_circuit:
        V_cap, V_gap = p2d.update_voltage_rc(dt, G_eff)

    f_current.write(f'{time:.6e} {J_tot:.6e} {J_displ:.6e} '
                    f'{V_gap:.6e} {V_cap:.6e}\n')
    f_current.flush()

    # Write output every N steps, with N = args.steps_per_output
    if step % args.steps_per_output == 0:
        i_output = step // args.steps_per_output
        p2d.write_solution(f'{args.siloname}_{i_output:04d}', i_output, time)

    streamers = [s for s in streamers if s.keep]
