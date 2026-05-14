# aseview

**Interactive molecular visualization for ASE (Atomic Simulation Environment)**

[![DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/kangmg/aseview)

aseview is a powerful, browser-based molecular viewer that integrates seamlessly with ASE. It provides interactive 3D visualization for molecules, crystals, and trajectories with support for normal mode animations.

## Features

- **Multiple Viewer Types**: Standard molecular viewer, overlay comparison, normal mode visualization, and interactive fragment selector
- **Rich Styling**: Multiple visual styles including Cartoon, Neon, Glossy, Metallic, Cinematic, Bubble, Grey, and 2D ([Styles](styles.md))
- **Themes**: Swap the entire look and feel with `dark`, `spring`, `glass`, `darkgreen`, `simple`, or custom themes ([Theming guide](theming.md))
- **Display Controls**: Hide hydrogens without deleting data, add render blur for saved images, and tune radius contrast
- **Trajectory Support**: Animate trajectories with synchronized energy and max-force plots
- **Normal Modes**: Visualize vibrational modes from ASE vibrations, ORCA Hessian files, or VASP OUTCAR files
- **Cell and Constraints**: Display periodic cells and highlight fixed atoms from `FixAtoms`
- **Radius Contrast**: Adjust element-radius differences from uniform radii to full element radii
- **Polyhedron**: Render coordination polyhedra for solid-state structures (CN ≥ 4)
- **Ring Highlight**: Fill aromatic and aliphatic rings with semi-transparent color
- **Fragment Selector**: Synchronized 2D/3D atom picker — select fragments visually, copy indices to clipboard
- **CLI & Python API**: Use from command line or Jupyter notebooks
- **Export**: Save as PNG, GIF, standalone HTML, or copy structures as `xyz`, `extxyz`, `cif`, and `POSCAR`

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
| **FragSelector** | Synchronized 2D/3D interactive atom selector | Fragment picking, QM/MM region assignment |

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
