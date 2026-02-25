# aseview

Molecular structure viewer for ASE (Atomic Simulation Environment).

**Status: pre-alpha** &nbsp;·&nbsp; [![DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/kangmg/aseview) &nbsp;·&nbsp; [📚 Documentation](https://kangmg.github.io/aseview)

## Installation

```bash
pip install aseview
```

For development:

```bash
git clone https://github.com/kangmg/aseview.git
cd aseview
pip install -e .
```

## CLI Usage

```bash
# Basic usage
aseview molecule.xyz

# View trajectory (all frames)
aseview trajectory.xyz -i :

# View specific frames
aseview trajectory.xyz -i 0:10    # frames 0-9
aseview trajectory.xyz -i -1      # last frame
aseview trajectory.xyz -i ::2     # every 2nd frame

# Specify file format
aseview POSCAR -f vasp

# Custom port
aseview molecule.xyz -p 9000

# Overlay multiple structures
aseview reactant.xyz product.xyz

# Overlay with colormap
aseview trajectory.xyz -v overlay --cmap viridis

# Normal mode visualization with ORCA Hessian
aseview molecule.xyz --hess orca.hess

# Save as HTML file
aseview molecule.xyz -o output.html

# Kill existing server on port
aseview molecule.xyz -k

# Help
aseview -h
```

## SSH Port Forwarding

When running on a remote server (e.g., HPC cluster, Docker container):

```bash
# 1. On remote server
aseview molecule.xyz -p 8080

# 2. On local machine (separate terminal)
ssh -L 8080:localhost:8080 user@remote-server

# 3. Open in local browser
# http://localhost:8080
```

For Docker with custom SSH port:

```bash
# Connect with port forwarding
ssh user@localhost -p 10011 -L 8080:localhost:8080

# Then run aseview inside container
aseview molecule.xyz -p 8080
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
    backgroundColor='#000000',
    showPolyhedron=True,    # coordination polyhedra (solid-state, CN >= 4)
    polyhedronOpacity=0.25,
    showRings=True,         # ring face highlighting (molecules, 4-8 atom rings)
    ringOpacity=0.35,
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

## JavaScript Module

Use aseview in any web page without Python:

```html
<div id="viewer" style="width:100%; height:500px;"></div>

<script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
<script>
    const viewer = new ASEView.MolecularViewer('#viewer');
    viewer.setData({
        symbols: ['O', 'H', 'H'],
        positions: [
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469]
        ]
    });
</script>
```

See the [JavaScript Module documentation](https://kangmg.github.io/aseview/js-module/) for full API reference.

## Supported Formats

All formats supported by ASE: xyz, cif, pdb, POSCAR, extxyz, etc.

## License

MIT
