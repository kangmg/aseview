# aseview

Molecular structure viewer for ASE (Atomic Simulation Environment).

[![PyPI](https://img.shields.io/pypi/v/aseview)](https://pypi.org/project/aseview/) &nbsp;·&nbsp; [![Docs](https://img.shields.io/badge/docs-online-blue)](https://kangmg.github.io/aseview) &nbsp;·&nbsp; [![Playground](https://img.shields.io/badge/playground-try_it-orange)](https://kangmg.github.io/aseview/playground/) &nbsp;·&nbsp; [![DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/kangmg/aseview)

## Features

- Interactive structure, trajectory, overlay, normal-mode, and fragment-selection viewers
- ASE-backed Python API and CLI, plus a browser-only JavaScript module
- Cell/PBC display, bonds, hydrogen bonds, charges, forces, and fixed-atom constraint highlighting
- Energy and max-force (Fmax) trajectory plots with independent force-plot toggling
- Radius Contrast control for reducing or restoring element-radius size differences
- Clipboard export for current frame or full trajectory as `xyz`, `extxyz`, `cif`, or `POSCAR`

## Installation

Recommended:

```bash
uv venv -p 3.11
uv pip install aseview
uv run aseview -h
```

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

# Normal mode visualization with ORCA Hessian or VASP OUTCAR
aseview molecule.xyz --hess orca.hess
aseview POSCAR --hess OUTCAR

# Save as HTML file
aseview molecule.xyz -o output.html

# Choose a visual theme
aseview molecule.xyz --theme spring
aseview molecule.xyz -t glass -o out.html

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

#### From VASP OUTCAR

```python
from ase.io import read
from aseview import NormalViewer

atoms = read('POSCAR')
viewer = NormalViewer.from_vasp(atoms, 'OUTCAR')
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
    style='neon',           # default, cartoon, neon, glossy, metallic, cinematic, rowan, grey
    bondThreshold=1.2,      # bond detection scale factor
    atomSize=0.5,
    radiusContrast=0.8,     # 0.0 = uniform radii, 1.0 = element radii
    radiusContrastMode='log',  # log or linear interpolation
    showCell=False,
    backgroundColor='#000000',
    showEnergyPlot=True,
    showForceMaxPlot=True,
    showConstraint=True,    # highlight FixAtoms constraints
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

### Fragment Selector (Interactive Atom Selection)

```python
from ase.build import molecule
from aseview import FragSelector

atoms = molecule('CH3OH')  # or read('molecule.xyz')
FragSelector(atoms).show()
```

The Fragment Selector shows a **synchronized 2D + 3D view** for picking atom subsets.
Selected atoms are highlighted in yellow in both panels.

Interaction modes (keyboard shortcuts: **N** / **C** / **R** / **L**):

| Mode | How to activate | 2D action |
|------|----------------|-----------|
| **Navigate** | N key | drag to pan, scroll to zoom |
| **Click** | C key | click individual atoms to toggle |
| **Rect** | R key | drag a rectangle to select atoms |
| **Lasso** | L key | draw a freeform region to select atoms |

Hold **Shift** while releasing to *add* to the existing selection instead of replacing it.

```python
# Save selected indices after interactive picking:
viewer = FragSelector(atoms, bondThreshold=1.2, style='cartoon')
viewer.save_html('selector.html')

# CLI usage (opens in browser):
# aseview molecule.xyz -v frag
```

## Viewer Types

| Viewer | Description |
|--------|-------------|
| MolecularViewer | Single structure or trajectory animation |
| NormalViewer | Normal mode vibration visualization |
| OverlayViewer | Compare multiple structures overlaid |
| FragSelector | Interactive 2D+3D atom selection with rect/lasso |

## VS Code Extension

The `vscode-extension/` package provides a native VS Code custom editor for structure files. It uses `ase-ts` for parsing and writing, so it does not require a Python runtime inside VS Code.

Supported default file patterns include `*.xyz`, `*.extxyz`, `*.cif`, `*.pdb`, `*.vasp`, `POSCAR`, and `CONTCAR`.

```bash
cd vscode-extension
npm install
npm run compile
npm run package
```

Install the generated `.vsix` from VS Code with **Extensions -> Install from VSIX...**.
On GitHub releases, the VS Code workflow builds the same package and attaches `aseview.vsix` to the release assets.

## Themes

aseview ships with multiple visual themes. Each theme is a complete HTML template set that controls the viewer's colour scheme, background, and UI style.

| Theme | Description |
|-------|-------------|
| `dark` | Default dark theme with deep grey background |
| `darkgreen` | Dark theme with green accent colours |
| `simple` | Minimal, low-distraction theme |
| `spring` | Light pastel theme with bright, airy colours |
| `glass` | Frosted-glass aesthetic with translucent UI panels |

### CLI

```bash
aseview molecule.xyz --theme spring
aseview molecule.xyz -t glass -o out.html
```

### Python API

```python
import aseview

# Per-viewer
viewer = aseview.MolecularViewer(atoms, theme='spring')
viewer = aseview.OverlayViewer([a, b], theme='glass')
viewer = aseview.NormalViewer(atoms, theme='dark')
viewer = aseview.FragSelector(atoms, theme='spring')

# Global — all subsequent viewers use this theme
aseview.set_theme('spring')
viewer1 = aseview.MolecularViewer(atoms)   # spring
viewer2 = aseview.FragSelector(atoms)       # spring

# Inspect
aseview.list_themes()   # ['dark', 'darkgreen', 'glass', 'simple', 'spring']
aseview.get_theme()     # current default
```

### JavaScript Module

```javascript
// Per-viewer
const v = new ASEView.MolecularViewer('#container', { theme: 'spring' });

// Global
ASEView.setTheme('glass');
const v2 = new ASEView.OverlayViewer('#container');  // glass
```

See the [Theming guide](https://kangmg.github.io/aseview/theming/) for instructions on creating custom themes.

## JavaScript Module

Use aseview in any web page without Python:

```html
<div id="viewer" style="width:100%; height:500px;"></div>

<script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
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

Input files are handled through ASE in Python and `ase-ts` in the VS Code extension. Common formats include `xyz`, `extxyz`, `cif`, `pdb`, `vasp`, `POSCAR`, and `CONTCAR`.

Viewer clipboard export supports:

| Export format | Notes |
|---------------|-------|
| `xyz` | Simple coordinates for current frame or trajectory |
| `extxyz` | Includes cell metadata and fixed/move-mask constraint data when present |
| `cif` | Includes cell/periodic structure data; constraints are ignored because CIF has no standard FixAtoms field |
| `POSCAR` | Includes cell and selective-dynamics flags for fixed/move-mask constraints |

## License

MIT
