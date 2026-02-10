# Quick Start

## Command Line Usage

### View a Single Structure

```bash
aseview2 molecule.xyz
```

This opens your default browser with an interactive 3D viewer.

### View a Trajectory

```bash
# View all frames
aseview2 trajectory.xyz

# View specific frames
aseview2 trajectory.xyz -i 0:10    # First 10 frames
aseview2 trajectory.xyz -i -1      # Last frame only
aseview2 trajectory.xyz -i ::2     # Every 2nd frame
```

### Compare Structures (Overlay)

```bash
# Overlay multiple files
aseview2 reactant.xyz product.xyz

# Overlay with colormap
aseview2 trajectory.xyz -v overlay --cmap viridis
```

### Normal Mode Visualization

```bash
# With ORCA Hessian file
aseview2 molecule.xyz --hess orca.hess
```

## Python API Usage

### In Jupyter Notebook

```python
from ase.io import read
from aseview import MolecularViewer

# Load structure
atoms = read("molecule.xyz")

# Create and display viewer
viewer = MolecularViewer(atoms)
viewer.show()
```

### Save as HTML

```python
viewer = MolecularViewer(atoms)
viewer.save_html("molecule.html")
```

### Customize Appearance

```python
viewer = MolecularViewer(
    atoms,
    style="neon",
    atomSize=0.5,
    bondThickness=0.15,
    backgroundColor="#000000"
)
viewer.show()
```

## SSH Remote Usage

When running on a remote server:

```bash
# On server
aseview2 molecule.xyz -p 8080

# On local machine
ssh -L 8080:localhost:8080 user@remote

# Open in local browser
# http://localhost:8080
```
