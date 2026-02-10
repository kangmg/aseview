# Overlay Comparison

Compare multiple molecular structures by overlaying them in a single view.

## Comparing Two Structures

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
    # Auto-detected as overlay when multiple files given
    aseview2 reactant.xyz product.xyz
    ```

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

### Color by Colormap (Gradient)

Gradient coloring for trajectories:

```python
viewer = OverlayViewer(
    trajectory,
    colorBy="Colormap",
    colormap="viridis"
)
```

Available colormaps:

| Colormap | Description |
|----------|-------------|
| `viridis` | Blue → Green → Yellow |
| `plasma` | Purple → Orange → Yellow |
| `coolwarm` | Blue → White → Red |
| `jet` | Rainbow spectrum |
| `rainbow` | Full color spectrum |
| `grayscale` | Black → White |

## Centering Molecules

By default, molecules are displayed at their original coordinates. To center each molecule at the origin for better comparison:

```python
viewer = OverlayViewer(
    structures,
    centerMolecules=True
)
```

This is useful when comparing conformers that are not pre-aligned.

## Interactive Controls

### Molecule Panel

For each overlaid molecule:

- **Visibility toggle**: Show/hide individual structures
- **Opacity slider**: Adjust transparency (0-100%)
- **Color picker**: Change molecule color

### Animation

When overlaying trajectory frames:

- **Play/Pause**: Animate through frames
- **Frame slider**: Jump to specific frame

## Examples

### Reaction Path

```python
from ase.io import read
from aseview import OverlayViewer

# Load IRC trajectory
irc = read("irc.xyz", index=":")

viewer = OverlayViewer(
    irc,
    colorBy="Colormap",
    colormap="coolwarm"  # Blue (reactant) → Red (product)
)
viewer.show()
```

### Conformer Comparison

```python
# Compare multiple conformers
conformers = [read(f"conf{i}.xyz") for i in range(1, 5)]

viewer = OverlayViewer(
    conformers,
    centerMolecules=True,  # Align at center
    colorBy="Molecule"
)
viewer.show()
```

### Optimization Overlay

```python
# Overlay first and last frames of optimization
traj = read("optimization.xyz", index=":")
first_last = [traj[0], traj[-1]]

viewer = OverlayViewer(
    first_last,
    colorBy="Molecule"
)
viewer.show()
```

### CLI with Colormap

```bash
# Overlay trajectory with plasma colormap
aseview2 trajectory.xyz -v overlay --cmap plasma

# Multiple files with viridis
aseview2 *.xyz --cmap viridis
```
