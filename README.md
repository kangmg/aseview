# aseview

A molecular viewer for ASE (Atomic Simulation Environment) data with React-based UI components.

## Features

- Interactive 3D molecular visualization
- Multiple viewing modes (Cartoon, Bubble, Glossy, Neon)
- Animation support
- Normal mode visualization
- Overlay visualization for comparing multiple molecules

## Installation

```bash
pip install aseview
```

## Usage

```python
import aseview
from ase.build import molecule

# Create a molecule
atoms = molecule('H2O')

# Visualize the molecule
viewer = aseview.MolecularViewer(atoms)
viewer.show()
```

## Development

To run the development version:

1. Clone the repository
2. Install dependencies:
   ```bash
   pip install -e .
   ```
3. Build the React components:
   ```bash
   cd ui
   npm install
   npm run build
   ```

## License

MIT