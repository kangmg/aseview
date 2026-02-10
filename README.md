# aseview

Molecular structure viewer for ASE (Atomic Simulation Environment).

**Status: pre-alpha**

## Installation

```bash
pip install -e .
```

## CLI Usage

```bash
# Basic usage
aseview2 molecule.xyz

# View trajectory (all frames)
aseview2 trajectory.xyz -i :

# View specific frames
aseview2 trajectory.xyz -i 0:10    # frames 0-9
aseview2 trajectory.xyz -i -1      # last frame
aseview2 trajectory.xyz -i ::2     # every 2nd frame

# Specify file format
aseview2 POSCAR -f vasp

# Custom port
aseview2 molecule.xyz -p 9000

# Overlay multiple structures
aseview2 reactant.xyz product.xyz

# Overlay with colormap
aseview2 trajectory.xyz -v overlay --cmap viridis

# Normal mode visualization with ORCA Hessian
aseview2 molecule.xyz --hess orca.hess

# Save as HTML file
aseview2 molecule.xyz -o output.html

# Kill existing server on port
aseview2 molecule.xyz -k

# Help
aseview2 -h
```

## SSH Port Forwarding

When running on a remote server (e.g., HPC cluster, Docker container):

```bash
# 1. On remote server
aseview2 molecule.xyz -p 8080

# 2. On local machine (separate terminal)
ssh -L 8080:localhost:8080 user@remote-server

# 3. Open in local browser
# http://localhost:8080
```

For Docker with custom SSH port:

```bash
# Connect with port forwarding
ssh user@localhost -p 10011 -L 8080:localhost:8080

# Then run aseview2 inside container
aseview2 molecule.xyz -p 8080
```

## Jupyter Notebook

### Quick Start

```python
from ase.io import read
from aseview import MolecularViewer

atoms = read('molecule.xyz')
viewer = MolecularViewer(atoms)
viewer.show()
```

### With Trajectory

```python
from ase.io import read
from aseview import MolecularViewer

# Read all frames
trajectory = read('trajectory.xyz', index=':')
viewer = MolecularViewer(trajectory)
viewer.show()
```

### Overlay Multiple Structures

```python
from ase.io import read
from aseview import OverlayViewer

reactant = read('reactant.xyz')
product = read('product.xyz')
viewer = OverlayViewer([reactant, product])
viewer.show()
```

### Overlay with Colormap

```python
from ase.io import read
from aseview import OverlayViewer

trajectory = read('optimization.xyz', index=':')
viewer = OverlayViewer(
    trajectory,
    colorBy='Colormap',
    colormap='viridis'  # viridis, plasma, coolwarm, jet, rainbow, grayscale
)
viewer.show()
```

### Align Molecules (RMSD Minimization)

```python
from ase.io import read
from aseview import OverlayViewer

structures = [read(f'conf{i}.xyz') for i in range(5)]
viewer = OverlayViewer(
    structures,
    alignMolecules=True,  # Kabsch rotation + Hungarian reordering
    colorBy='Molecule'
)
viewer.show()
```

### Normal Mode Visualization

#### From ASE Vibrations

```python
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from aseview import NormalViewer

# Create or load molecule
atoms = Atoms('H2O', positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
atoms.calc = EMT()

# Optimize structure
opt = BFGS(atoms)
opt.run(fmax=0.01)

# Calculate vibrations
vib = Vibrations(atoms, name='vib')
vib.run()
vib.summary()

# Visualize normal modes
viewer = NormalViewer(atoms, vibrations=vib)
viewer.show()
```

#### From ORCA Hessian File

```python
from ase.io import read
from aseview import NormalViewer

atoms = read('molecule.xyz')
viewer = NormalViewer.from_orca(atoms, 'orca.hess')
viewer.show()
```

Features:
- Mode selector dropdown with frequencies
- Sinusoidal animation of atomic displacements
- Amplitude slider to control displacement magnitude
- Show Vectors toggle to display mode displacement arrows
- Imaginary frequencies (transition states) shown in red

### Using view_molecule Helper

```python
from aseview.jupyter import view_molecule
from ase.io import read

atoms = read('molecule.xyz')
view_molecule(atoms, viewer_type='molecular', height=600)
```

### Custom Settings

```python
from aseview import MolecularViewer

viewer = MolecularViewer(
    atoms,
    style='neon',           # default, cartoon, neon, glossy, metallic, rowan, grey
    bondThreshold=1.2,      # bond detection scale factor
    atomSize=0.5,
    showCell=False,
    backgroundColor='#000000'
)
viewer.show(width='100%', height=800)
```

### Save to HTML

```python
viewer.save_html('output.html')
```

## Viewer Types

| Viewer | Description |
|--------|-------------|
| MolecularViewer | Single structure or trajectory animation |
| NormalViewer | Normal mode vibration visualization |
| OverlayViewer | Compare multiple structures overlaid |

## Supported Formats

All formats supported by ASE: xyz, cif, pdb, POSCAR, extxyz, etc.

## License

MIT
