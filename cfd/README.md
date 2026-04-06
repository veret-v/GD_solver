# CFD Solver

2D Euler solver on unstructured triangular meshes, organized to match the HTML lecture structure.

## Dependency

The project uses the official Gmsh C++ API to load `.msh` files.

Set `GMSH_DIR` to the root of the Gmsh SDK:

- Windows example: `C:\sdk\gmsh`
- Linux example: `/opt/gmsh-sdk`

Expected layout:

- `${GMSH_DIR}/include/gmsh.h`
- `${GMSH_DIR}/lib/gmsh.lib` or `${GMSH_DIR}/lib/libgmsh.so`

## Build

Windows:

```powershell
cd cfd
cmake -S . -B build -DGMSH_DIR="C:/path/to/gmsh-sdk"
cmake --build build --config Release
```

Linux/macOS:

```bash
cd cfd
cmake -S . -B build -DGMSH_DIR=/path/to/gmsh-sdk
cmake --build build -j
```

## Example: Cylinder In Channel

Generate the mesh:

```bash
gmsh -2 meshes/cylinder_channel.geo -o meshes/cylinder_channel.msh
```

Run the solver:

```bash
build/cfd --mesh-type gmsh --mesh meshes/cylinder_channel.msh --mach 0.3 --cfl 0.4 --tend 5.0 --io 50 --time rk2 --scheme godunov --cyl-x 1.5 --cyl-y 1.0 --cyl-r 0.25
```
