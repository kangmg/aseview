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

<iframe src="../assets/viewers/water.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

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

<iframe src="../assets/viewers/ethanol_neon.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Benzene (Glossy Style)

```python
from ase.build import molecule
from aseview import MolecularViewer

benzene = molecule('C6H6')
viewer = MolecularViewer(benzene, style="glossy")
viewer.show()
```

<iframe src="../assets/viewers/benzene_glossy.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Caffeine (Metallic Style)

```python
from aseview import MolecularViewer

viewer = MolecularViewer(caffeine, style="metallic")
viewer.show()
```

<iframe src="../assets/viewers/caffeine.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

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
