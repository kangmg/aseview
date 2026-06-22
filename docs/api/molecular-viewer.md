# MolecularViewer

Standard 3D molecular viewer with support for trajectories and animations.

## Constructor

```python
MolecularViewer(data, **kwargs)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `data` | `Atoms`, `List[Atoms]`, `Dict`, `str` | Molecular data to visualize | Required |
| `theme` | `str` | Visual theme (`"dark"`, `"spring"`, `"glass"`, `"darkgreen"`, `"simple"`, …). `None` uses the global default. | `None` |

### Keyword Arguments (Settings)

#### Geometry Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atomSize` | `float` | Atom sphere radius scale | `0.4` |
| `bondThickness` | `float` | Bond cylinder radius | `0.09` |
| `bondThreshold` | `float` | Bond detection threshold (multiplier for covalent radii sum) | `1.2` |
| `radiusContrast` | `float` | Relative atom radius contrast (`0.0` = uniform radii, `1.0` = full element radii) | `1.0` |
| `radiusContrastMode` | `str` | Radius contrast mapping: `"linear"` or `"log"` | `"log"` |

#### Style Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `style` | `str` | Visual style name | `"Cartoon"` |
| `backgroundColor` | `str` | Background color (hex) | `"#1f2937"` |

##### Available Styles

| Value | Description |
|-------|-------------|
| `"default"` | Standard CPK coloring |
| `"2d"` | Flat, high-contrast 2D-style rendering |
| `"cartoon"` | Cartoon style with black bonds |
| `"neon"` | Glowing neon effect |
| `"glossy"` | Shiny reflective surface |
| `"metallic"` | Metallic appearance |
| `"cinematic"` | Cinematic glossy ball-and-stick style |
| `"rowan"` | Rowan-inspired style |
| `"grey"` | Greyscale rendering with bonds matching the active atom colors |
| `"bubble"` | Bubble-like appearance |

See [Styles](../styles.md) for interactive previews of every style.

#### Display Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `showCell` | `bool` | Show unit cell | `True` |
| `showBond` | `bool` | Show covalent bonds | `True` |
| `hideHydrogens` | `bool` | Hide hydrogen atoms and their bonds without removing them from the data | `False` |
| `showBlur` | `bool` | Apply blur to the rendered molecule canvas | `False` |
| `blurStrength` | `float` | Blur radius in pixels when `showBlur` is enabled | `1.5` |
| `showHBond` | `bool` | Show hydrogen bonds | `False` |
| `hBondThreshold` | `float` | H···A distance cutoff (Å) for hydrogen-bond detection | `2.5` |
| `showShadow` | `bool` | Enable shadows | `False` |
| `showEnergyPlot` | `bool` | Show energy plot (if energy data available) | `False` |
| `showForceMaxPlot` | `bool` | Include max force (`fmax`) in the trajectory plot when forces are available | `True` |
| `showForces` | `bool` | Show force vectors | `False` |
| `showShading` | `bool` | Enable 3D shading effect | `True` |

#### Color Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `colorBy` | `str` | Atom coloring mode: `"Element"`, `"Charge"`, or `"Magmom"` | `"Element"` |
| `colorScheme` | `str` | Element color scheme: `"Jmol"`, `"CPK"`, `"PyMOL"`, or `"VMD"` | `"Jmol"` |
| `normalizeCharges` | `bool` | Normalize charges to symmetric range | `False` |
| `chargeColormap` | `str` | Colormap for charge visualization | `"coolwarm"` |
| `showConstraint` | `bool` | Highlight fixed atoms with a yellow overlay | `False` |

!!! note "Charge Coloring"
    The `colorBy="Charge"` option is only available when atoms have charge data.
    Charges can be set via:

    - `atoms.arrays['charges']` - per-atom array
    - `atoms.info['charges']` - list in info dict

    The colormap uses a diverging scheme:

    - **Blue** → positive charges
    - **White** → neutral (zero)
    - **Red** → negative charges

!!! note "Magmom Coloring"
    The `colorBy="Magmom"` option is only available when atoms have per-atom magnetic moment data.
    Magnetic moments can be set via:

    - `atoms.arrays['magmoms']` or `atoms.arrays['initial_magmoms']` - per-atom arrays
    - `atoms.info['magmom']` or `atoms.info['magmoms']` - lists in the info dict

    Magmom coloring does not use the `normalizeCharges` setting. Values are shown on a zero-centered diverging scale, with color intensity based on distance from zero.

!!! note "Constraint Visualization"
    The `showConstraint=True` option is only available when atoms have `FixAtoms`
    constraints. Fixed atoms are rendered with a semi-transparent yellow sphere overlay
    on top of their normal element color.

    ```python
    from ase.constraints import FixAtoms

    atoms.set_constraint(FixAtoms(indices=[0, 1, 2]))
    viewer = MolecularViewer(atoms, showConstraint=True)
    ```

!!! note "Radius Contrast"
    `radiusContrast` controls how different the element radii appear relative to one another.
    Use `0.0` for uniform atom sizes and `1.0` for the built-in element radii.
    `radiusContrastMode="log"` compresses large radius differences more strongly than `"linear"`.

    ```python
    viewer = MolecularViewer(atoms, radiusContrast=0.6, radiusContrastMode="log")
    ```

!!! note "Trajectory Plot"
    `showEnergyPlot=True` displays energy when each frame provides `atoms.get_potential_energy()`
    or `atoms.info['energy']`. If forces are present, the plot can also show per-frame
    max force (`fmax`) with `showForceMaxPlot=True`.

!!! note "Hydrogen Bonds"
    `showHBond=True` detects and renders hydrogen bonds as dashed lines.
    Criteria: donor H bonded to N/O/F, acceptor N/O/F, H···A ≤ `hBondThreshold` (default 2.5 Å),
    D···A ≤ 3.5 Å, with direct-bonded (1,2) and shared-neighbor (1,3) donor–acceptor pairs excluded.

#### Force Vector Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `forceScale` | `float` | Scale factor for force arrows | `1.0` |

#### Animation Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `animationSpeed` | `int` | Animation speed (ms per frame) | `30` |

## Clipboard Copy Formats

The animation toolbar can copy the current frame or the full loaded trajectory.
The copy-format chip cycles through `xyz`, `extxyz`, `cif`, and `POSCAR`.

| Format | Cell Support | Constraint Support |
|--------|--------------|--------------------|
| `xyz` | No | No |
| `extxyz` | Yes | Fixed/move-mask constraint metadata |
| `cif` | Yes | No — CIF has no standard fixed-atom field |
| `POSCAR` | Yes | Selective-dynamics flags for fixed/move-mask constraints |

#### View Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `viewMode` | `str` | Camera projection: `"Perspective"` or `"Orthographic"` | `"Perspective"` |
| `viewPreset` | `str`, `None` | Initial named camera direction (`"top-c"`, `"side-a"`, `"front"`, etc.) | `None` |
| `viewDirection` | `list[float]`, `None` | Explicit target-to-camera direction vector | `None` |
| `viewEuler` | `list[float]`, `None` | `[rx, ry, rz]` degrees in XYZ order, applied to `[0, 0, 1]` | `None` |
| `viewUp` | `list[float]`, `None` | Optional camera up-vector hint | `None` |
| `viewFit` | `float` | Camera fit multiplier | `1.0` |
| `rotationMode` | `str` | Rotation control: `"TrackBall"` or `"Orbit"` | `"TrackBall"` |
| `selectionMode` | `str` | Selection mode: `"Lasso"` or `"Click"` | `"Lasso"` |

`viewDirection` and preset directions are target-to-camera vectors. Presets
`top`, `bottom`, `front`, `back`, `left`, and `right` use Cartesian axes.
Presets `top-c`, `bottom-c`, `side-a`, and `side-b` use valid unit-cell
vectors. Aliases `c` and `top` prefer the cell `c` axis, while `a` and `b`
prefer the cell `a` and `b` axes. Missing, non-periodic, zero-length, or
degenerate cells fall back to Cartesian directions.

## Methods

### show()

Display the viewer in a Jupyter notebook.

```python
viewer.show(width='100%', height=600)
```

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `width` | `str`, `int` | Width of viewer | `'100%'` |
| `height` | `int` | Height in pixels | `600` |

### get_html()

Get the HTML content as a string.

```python
html = viewer.get_html()
```

**Returns:** `str` - Complete HTML document

### save_html()

Save the viewer as an HTML file.

```python
viewer.save_html(filename)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `filename` | `str` | Output file path |

!!! note "Image export"
    `MolecularViewer` can save PNG/GIF images from Python when the optional
    export dependency is installed:

    ```bash
    pip install "aseview[export]"
    python -m playwright install chromium
    ```

    ```python
    viewer.save_png(
        "molecule.png",
        width=1600,
        height=1200,
        transparent=False,
        background_color="#ffffff",
    )
    viewer.save_gif("trajectory.gif", frames=30, quality="high")
    ```

    PNG quality is controlled by `scale`, `width`, and `height`; GIF
    `quality="low" | "medium" | "high"` maps to encoder sampling quality.
    The CLI does not provide `--save-png` or `--save-gif`.

    In the browser JavaScript module, `ASEView.MolecularViewer` exposes
    `setView(viewSpec)`, `resetView()`, `savePNG(options)`, and
    `saveGIF(options)`. Export options are `filename`, `download`,
    `returnDataUrl`, `scale`, `width`, `height`, `transparent`,
    `backgroundColor`, plus GIF-only `frames`, `delay`, and `sampleInterval`.
    Export promises resolve with
    `{ ok: true, type: "png" | "gif", filename, dataUrl?, width?, height? }` or
    reject with an `Error` carrying `code`, `message`, `type`, and when
    available `requestId`.

## Examples

### Basic Usage

```python
from ase.io import read
from aseview import MolecularViewer

atoms = read("molecule.xyz")
viewer = MolecularViewer(atoms)
viewer.show()
```

### Custom Styling

```python
viewer = MolecularViewer(
    atoms,
    style="neon",
    atomSize=0.5,
    bondThickness=0.15,
    backgroundColor="#000000",
    viewPreset="top-c",
    viewFit=1.1
)
viewer.show()
```

### Trajectory with Forces

```python
from ase.io import read

trajectory = read("md.xyz", index=":")
viewer = MolecularViewer(
    trajectory,
    showForces=True,
    forceScale=2.0,
    showEnergyPlot=True,
    showForceMaxPlot=True
)
viewer.show()
```

### Radius Contrast

```python
viewer = MolecularViewer(
    atoms,
    atomSize=0.45,
    radiusContrast=0.5,
    radiusContrastMode="linear",
)
viewer.show()
```

### Partial Charge Visualization

```python
import numpy as np
from ase.build import molecule
from aseview import MolecularViewer

# Create molecule with partial charges
water = molecule("H2O")

# Set charges (e.g., from DFT Mulliken/Hirshfeld analysis)
# O is negative, H atoms are positive
charges = np.array([-0.82, 0.41, 0.41])
water.arrays['charges'] = charges

# Visualize with charge coloring
viewer = MolecularViewer(
    water,
    colorBy="Charge",
    normalizeCharges=False  # Use actual charge values
)
viewer.show()
```

!!! tip "Setting Charges"
    Charges can also be set via `atoms.info['charges']`:
    ```python
    water.info['charges'] = [-0.82, 0.41, 0.41]
    ```

!!! tip "Setting Magnetic Moments"
    Magnetic moments can be set via `atoms.info['magmom']`:
    ```python
    water.info['magmom'] = [1.0, -1.0, 0.0]
    viewer = MolecularViewer(water, colorBy="Magmom")
    ```

### Constraint Visualization

```python
from ase.io import read
from ase.constraints import FixAtoms
from aseview import MolecularViewer

atoms = read("surface.xyz")
atoms.set_constraint(FixAtoms(indices=list(range(12))))  # fix bottom layer

viewer = MolecularViewer(atoms, showConstraint=True)
viewer.show()
```

### Hydrogen Bond Visualization

```python
from ase.io import read
from aseview import MolecularViewer

atoms = read("water_dimer.xyz")
viewer = MolecularViewer(atoms, showHBond=True, hBondThreshold=2.8)
viewer.show()
```

### Color Scheme

```python
# CPK, PyMOL, and VMD element color tables are available
viewer = MolecularViewer(atoms, colorScheme="PyMOL")
viewer.show()
```

### Save to File

```python
viewer = MolecularViewer(atoms, style="glossy")
viewer.save_html("molecule_viewer.html")
```
