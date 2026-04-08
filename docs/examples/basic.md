# Basic Visualization

## Simple Usage

=== "Python"

    ```python
    from ase.io import read
    from aseview import MolecularViewer

    atoms = read("molecule.xyz")
    viewer = MolecularViewer(atoms, style="cartoon")
    viewer.show()
    ```

=== "CLI"

    ```bash
    aseview molecule.xyz --style cartoon
    ```

<iframe src="../../assets/viewers/water.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Available Styles

| Style | Description |
|-------|-------------|
| `cartoon` | Cartoon style with black bonds (default) |
| `default` | Standard CPK coloring |
| `neon` | Glowing neon effect |
| `glossy` | Shiny reflective surface |
| `metallic` | Metallic appearance |
| `rowan` | Rowan-inspired style |
| `bubble` | Bubble-like appearance |
| `grey` | Greyscale rendering |

### Neon Style Example

```python
viewer = MolecularViewer(atoms, style="neon", backgroundColor="#000000")
viewer.show()
```

<iframe src="../../assets/viewers/ethanol_neon.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Partial Charge Visualization

Visualize partial charges using `colorBy="Charge"`. Colors: red (negative) → white (neutral) → blue (positive).

```python
import numpy as np
from ase.io import read
from aseview import MolecularViewer

atoms = read("molecule.xyz")
atoms.arrays['charges'] = np.array([...])  # partial charges

viewer = MolecularViewer(atoms, colorBy="Charge")
viewer.show()
```

<iframe src="../../assets/viewers/phosphine_charges.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

### Charge Options

| Option | Description |
|--------|-------------|
| `colorBy="Charge"` | Enable charge coloring |
| `normalizeCharges=True` | Scale colors to full range |
| `showChargeLabels=True` | Display charge values on atoms |

---

## Constraint Visualization

Highlight fixed atoms with `colorBy="Constraint"`. Fixed atoms (from `FixAtoms`) appear in **grey**; free atoms keep their element color. The option is only available when constraints are present.

```python
from ase.constraints import FixAtoms
from aseview import MolecularViewer

atoms.set_constraint(FixAtoms(indices=[0, 1, 2]))
viewer = MolecularViewer(atoms, colorBy="Constraint")
viewer.show()
```

See the [Solid-State examples](../examples/solid-state.md#live-demo-cu111-surface--fixed-atoms-fixatoms-constraint) for a full surface slab demo.

---

## Customization

```python
viewer = MolecularViewer(
    atoms,
    style="metallic",
    atomSize=0.5,           # Atom radius scale
    bondThickness=0.15,     # Bond radius
    bondThreshold=1.2,      # Bond detection threshold
    backgroundColor="#1a1a2e"
)
viewer.show()
```

## Saving to HTML

```python
viewer = MolecularViewer(atoms)
viewer.save_html("molecule.html")
```

Creates a standalone HTML file viewable in any browser without Python.
