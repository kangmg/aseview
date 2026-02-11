# Basic Visualization

This example demonstrates basic molecular visualization with different styles.

## Water Molecule (Cartoon Style)

=== "Python"

    ```python
    from ase.build import molecule
    from aseview import MolecularViewer

    water = molecule('H2O')
    viewer = MolecularViewer(water, style="cartoon")
    viewer.show()
    ```

=== "CLI"

    ```bash
    aseview2 water.xyz --style cartoon
    ```

<iframe src="../../assets/viewers/water.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Ethanol (Neon Style)

```python
from ase.build import molecule
from aseview import MolecularViewer

ethanol = molecule('CH3CH2OH')
viewer = MolecularViewer(
    ethanol,
    style="neon",
    backgroundColor="#000000"
)
viewer.show()
```

<iframe src="../../assets/viewers/ethanol_neon.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Benzene (Glossy Style)

```python
from ase.build import molecule
from aseview import MolecularViewer

benzene = molecule('C6H6')
viewer = MolecularViewer(benzene, style="glossy")
viewer.show()
```

<iframe src="../../assets/viewers/benzene_glossy.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Caffeine (Metallic Style)

```python
from aseview import MolecularViewer

viewer = MolecularViewer(caffeine, style="metallic")
viewer.show()
```

<iframe src="../../assets/viewers/caffeine.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Partial Charge Visualization

You can visualize partial charges using `colorBy="Charge"`. Atoms are colored using a coolwarm colormap: blue (negative) → white (neutral) → red (positive).

### Phosphine with Partial Charges

```python
import numpy as np
from ase import Atoms
from aseview import MolecularViewer

# Create phosphine molecule with partial charges
phosphine = Atoms(symbols=[...], positions=[...])  # 47 atoms
charges = [0.0, 0.0, 0.065, ...]  # partial charges from DFT
phosphine.arrays['charges'] = np.array(charges)

viewer = MolecularViewer(phosphine, style="cartoon", colorBy="Charge")
viewer.show()
```

<iframe src="../../assets/viewers/phosphine_charges.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

### Ethanol with Partial Charges

```python
from ase.build import molecule
import numpy as np
from aseview import MolecularViewer

ethanol = molecule("CH3CH2OH")

# Partial charges: O is negative, acidic H is positive
charges = [-0.18, 0.06, 0.06, 0.06, 0.15, 0.04, 0.04, -0.65, 0.42]
ethanol.arrays['charges'] = np.array(charges)

viewer = MolecularViewer(ethanol, style="glossy", colorBy="Charge")
viewer.show()
```

<iframe src="../../assets/viewers/ethanol_charges.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Available Styles

| Style | Description |
|-------|-------------|
| `default` | Standard CPK coloring |
| `cartoon` | Cartoon style with black bonds (default) |
| `neon` | Glowing neon effect |
| `glossy` | Shiny reflective surface |
| `metallic` | Metallic appearance |
| `rowan` | Rowan-inspired style |
| `grey` | Greyscale rendering |
| `bubble` | Bubble-like appearance |

## Customizing Appearance

```python
viewer = MolecularViewer(
    atoms,
    atomSize=0.5,           # Larger atoms
    bondThickness=0.15,     # Thicker bonds
    bondThreshold=1.2,      # More bonds detected
    backgroundColor="#000000",  # Black background
    style="metallic"
)
viewer.show()
```

## Saving to HTML

```python
viewer = MolecularViewer(atoms)
viewer.save_html("molecule.html")
```

This creates a standalone HTML file that can be opened in any browser without requiring Python.
