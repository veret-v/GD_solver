## Sod FLIC code bundle

This folder contains the minimal source set used for the Sod shock tube run with FLIC.

Files:
- `main.cpp` - program entry point and time stepping loop; FLIC dispatch is under `equation_type == 8`
- `flic.cpp`, `flic.h` - FLIC implementation
- `solver.cpp`, `solver.h` - CFL, state conversion, boundary handling, analytic solution helpers
- `grid.h` - shared variable indices and constants
- `sod_flic.ini` - run configuration used for the generated results

Markers required by the checklist:
- `CHECK: FLIC_LAGRANGE` in `flic.cpp`
- `CHECK: FLIC_EULER` in `flic.cpp`
- `CHECK: FLIC_CONSERV` in `flic.cpp`
- `CHECK: FLIC_CFL` in `solver.cpp`

Generated plots and values are stored in:
- `../01_sod_flic/`
