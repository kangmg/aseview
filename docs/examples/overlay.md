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
    aseview reactant.xyz product.xyz
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

## Atom Selection

Control which atoms are included in the overlay using `index_list`.
The mode is inferred automatically from the value type — no separate mode flag needed.

### all — use every atom (default)

```python
viewer = OverlayViewer([mol1, mol2, mol3])  # index_list=None
```

### same — identical indices for every structure

Useful when all structures share the same atom ordering and you want to focus on a subset:

```python
# Show only the first 3 heavy atoms across all conformers
viewer = OverlayViewer(conformers, index_list=[0, 1, 2])
viewer.show()
```

### indices — different indices per structure

Useful when atom ordering differs between structures (e.g. reactant vs. product with different numbering):

```python
# reactant: atoms 0–5 are the core; product: atoms 2–7 are the core
viewer = OverlayViewer(
    [reactant, product],
    index_list=[[0, 1, 2, 3, 4, 5], [2, 3, 4, 5, 6, 7]]
)
viewer.show()
```

!!! note
    In `indices` mode the length of `index_list` must equal the number of structures,
    otherwise a `ValueError` is raised.

---

## Initial Structure Visibility

For readability, overlays with many structures start with only the first three
structures visible. The molecule panel still lists every structure and lets you
toggle each one interactively.

```python
# Show all trajectory frames immediately
viewer = OverlayViewer(trajectory, all_visible=True)
viewer.show()
```

```python
# Show only first and last frames initially
viewer = OverlayViewer(
    trajectory,
    visible_indices=[0, len(trajectory) - 1],
)
viewer.show()
```

---

## Interactive Controls

| Control | Description |
|---------|-------------|
| **Visibility toggle** | Show/hide individual structures |
| **Opacity slider** | Adjust transparency |

---

## CLI Options

```bash
# With colormap
aseview trajectory.xyz -v overlay --cmap viridis

# Multiple files
aseview *.xyz --cmap plasma
```
