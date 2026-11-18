# Streamer velocity from 3D Silo frames 0000–0051

Input: `output/LongElectrode_3D__NNNN.silo`.

The saved Silo grids do not contain explicit streamer-head identities,
coordinates, or velocity vectors. The output estimates the leading conducting
front's axial velocity, not individual head speed or arc-length speed.

For each frame, find the greatest z coordinate of a real (non-ghost) cell center
where `sigma_tot` is at least 1e-3 S/m. Repeat at 1e-4 and 1e-2 S/m to assess
threshold sensitivity. The cell-center front is quantized at the local grid
spacing. Branch switching can influence this leading-front measurement.

Read physical times directly from the Silo root `dtime` variable. Calculate
`v_z = (z_front[i+1] - z_front[i]) / (time[i+1] - time[i])` and plot at the
interval midpoint. No smoothing or model velocity law is used. Frame zero has
no conducting channel, so its front and the first interval velocity are NaN.

Build and reproduce from the cocydimo_mychanges directory:

```bash
cc -O2 -Iafivo/external_libraries/silo/include analysis/streamer_velocity/extract.c afivo/external_libraries/silo/lib/libsilo.a -lm -o analysis/streamer_velocity/extract
analysis/streamer_velocity/extract output > analysis/streamer_velocity/fronts.csv
cc -O2 -Iafivo/external_libraries/silo/include analysis/streamer_velocity/read_times.c afivo/external_libraries/silo/lib/libsilo.a -lm -o analysis/streamer_velocity/read_times
analysis/streamer_velocity/read_times output > analysis/streamer_velocity/times.csv
/usr/bin/python analysis/streamer_velocity/plot_velocity.py
```

Files: PNG/PDF plot; head_positions.csv; head_velocity.csv; extraction sources.
