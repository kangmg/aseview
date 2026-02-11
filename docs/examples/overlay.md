# Overlay Comparison

Compare multiple molecular structures by overlaying them in a single view.

## Live Demo

<iframe src="../../assets/viewers/overlay_conformers.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Basic Usage

=== "Python"

    ```python
    from ase.io import read
    from aseview import OverlayViewer

    reactant = read("reactant.xyz")
    product = read("product.xyz")

    viewer = OverlayViewer([reactant, product])
    viewer.show()
    ```

=== "CLI"

    ```bash
    aseview2 reactant.xyz product.xyz
    ```

---

## Color Modes

### Color by Atom (Default)
Each atom colored by element (CPK colors):
```python
viewer = OverlayViewer(structures, colorBy="Atom")
```

### Color by Molecule
Each structure gets a distinct color:
```python
viewer = OverlayViewer(structures, colorBy="Molecule")
```

### Color by Colormap
Gradient coloring for trajectories:
```python
viewer = OverlayViewer(trajectory, colorBy="Colormap", colormap="viridis")
```

<iframe src="../../assets/viewers/overlay_colormap.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

### Available Colormaps

| Colormap | Description |
|----------|-------------|
| `viridis` | Blue → Green → Yellow |
| `plasma` | Purple → Orange → Yellow |
| `coolwarm` | Blue → White → Red |
| `jet` | Blue → Cyan → Yellow → Red |

---

## Molecule Alignment (Kabsch)

Align structures for optimal RMSD comparison:

```python
viewer = OverlayViewer(structures, alignMolecules=True)
```

This performs:
1. **Centering**: Both molecules centered at origin
2. **Hungarian reordering**: Match atoms by element and position
3. **Kabsch rotation**: Optimal rotation via SVD to minimize RMSD

---

## Interactive Controls

| Control | Description |
|---------|-------------|
| **Visibility toggle** | Show/hide individual structures |
| **Opacity slider** | Adjust transparency |
| **Align Molecules** | Toggle Kabsch alignment |

---

## CLI Options

```bash
# With colormap
aseview2 trajectory.xyz -v overlay --cmap viridis

# Multiple files
aseview2 *.xyz --cmap plasma
```
