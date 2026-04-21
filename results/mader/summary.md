# Checklist Bundle

Source outputs: /Users/sergei/CLionProjects/GD_solver/tests/mader_2de/output_audit

## sod
Sod checklist report
[PASS] grid: Nx=100, Ny=3
[PASS] gamma: gamma=1.4
[PASS] viscosity: mader_visc=0.5
[PASS] cfl: cfl=0.5
[PASS] wave positions: contact_x=0.7000, shock_x=0.8400
[FAIL] contact smearing: contact_width=12
[FAIL] shock smearing: shock_width=6
[FAIL] control value x≈0.60: rho=0.3851, p=0.3111, u=0.9889
[FAIL] control value x≈0.75: rho=0.3184, p=0.3117, u=1.0044

## cj
CJ checklist report
[PASS] grid: Nx=200, Ny=3
[PASS] gamma: gamma=3.0
[PASS] viscosity: mader_visc=0.5
[PASS] cfl: cfl=0.5
[FAIL] arrhenius rate: reaction_rate=500000000.0
[PASS] min temperature: min_temperature=1000.0
[PASS] gasw threshold: gasw_threshold=0.02
[PASS] front speed: D=0.8919
[PASS] front linearity: r2=0.9976
[FAIL] reaction zone width: zone_width=1
[FAIL] pressure plateau: mean=0.1861, std=0.0395

## corner
Corner checklist report
[PASS] domain: x_max=7.0, y_max=4.0
[PASS] grid: Nx=140, Ny=80
[PASS] channel geometry: step_x_end=2.0, step_y_end=3.0
[PASS] in-channel speed: D=0.8703, r2=0.9937
[PASS] curved front after exit: curvature=2.4500
[PASS] dead zone contrast: p_top=0.0768, p_bottom=0.0000, w_top=0.8341, w_bottom=1.0000
[PASS] stability: rho_min=1.156180, p_min=3.680000e-08
[PASS] shargatov implementation present, runtime flag: shargatov_correction=0

