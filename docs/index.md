# aseview

**Interactive molecular visualization for ASE (Atomic Simulation Environment)**

aseview is a powerful, browser-based molecular viewer that integrates seamlessly with ASE. It provides interactive 3D visualization for molecules, crystals, and trajectories with support for normal mode animations.

## Features

- **Multiple Viewer Types**: Standard molecular viewer, overlay comparison, and normal mode visualization
- **Rich Styling**: Multiple visual styles including Cartoon, Neon, Glossy, Metallic, and more
- **Trajectory Support**: Animate molecular dynamics trajectories with energy plot synchronization
- **Normal Modes**: Visualize vibrational modes from ORCA Hessian files
- **CLI & Python API**: Use from command line or Jupyter notebooks
- **Export**: Save as PNG, GIF, or standalone HTML files

## Quick Example

=== "Python"

    ```python
    from ase.io import read
    from aseview import MolecularViewer

    atoms = read("molecule.xyz")
    viewer = MolecularViewer(atoms)
    viewer.show()  # In Jupyter notebook
    ```

=== "CLI"

    ```bash
    aseview molecule.xyz
    ```

## Viewer Types

| Viewer | Description | Use Case |
|--------|-------------|----------|
| **MolecularViewer** | Standard 3D molecular viewer with animation | Single structures, trajectories |
| **OverlayViewer** | Overlay multiple structures for comparison | Comparing conformers, reaction paths |
| **NormalViewer** | Vibrational mode visualization | Frequency analysis, IR/Raman |

## Installation

```bash
pip install aseview
```

Or install from source:

```bash
git clone https://github.com/kangmg/aseview.git
cd aseview
pip install -e .
```
