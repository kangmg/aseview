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

### Keyword Arguments (Settings)

#### Geometry Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atomSize` | `float` | Atom sphere radius scale | `0.4` |
| `bondThickness` | `float` | Bond cylinder radius | `0.1` |
| `bondThreshold` | `float` | Bond detection threshold (multiplier for covalent radii sum) | `1.0` |

#### Style Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `style` | `str` | Visual style name | `"Cartoon"` |
| `backgroundColor` | `str` | Background color (hex) | `"#1f2937"` |

##### Available Styles

| Value | Description |
|-------|-------------|
| `"default"` | Standard CPK coloring |
| `"cartoon"` | Cartoon style with black bonds |
| `"neon"` | Glowing neon effect |
| `"glossy"` | Shiny reflective surface |
| `"metallic"` | Metallic appearance |
| `"rowan"` | Rowan-inspired style |
| `"grey"` | Greyscale rendering |
| `"bubble"` | Bubble-like appearance |

#### Display Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `showCell` | `bool` | Show unit cell | `True` |
| `showBond` | `bool` | Show bonds | `True` |
| `showShadow` | `bool` | Enable shadows | `False` |
| `showEnergyPlot` | `bool` | Show energy plot (if energy data available) | `False` |
| `showForces` | `bool` | Show force vectors | `False` |
| `showShading` | `bool` | Enable 3D shading effect | `True` |

#### Color Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `colorBy` | `str` | Atom coloring mode: `"Element"` or `"Charge"` | `"Element"` |
| `normalizeCharges` | `bool` | Normalize charges to symmetric range | `False` |
| `chargeColormap` | `str` | Colormap for charge visualization | `"coolwarm"` |

!!! note "Charge Coloring"
    The `colorBy="Charge"` option is only available when atoms have charge data.
    Charges can be set via:

    - `atoms.arrays['charges']` - per-atom array
    - `atoms.info['charges']` - list in info dict

    The colormap uses a diverging scheme:

    - **Blue** → negative charges
    - **White** → neutral (zero)
    - **Red** → positive charges

#### Force Vector Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `forceScale` | `float` | Scale factor for force arrows | `1.0` |

#### Animation Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `animationSpeed` | `int` | Animation speed (ms per frame) | `30` |

#### View Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `viewMode` | `str` | Camera projection: `"Perspective"` or `"Orthographic"` | `"Perspective"` |
| `rotationMode` | `str` | Rotation control: `"TrackBall"` or `"Orbit"` | `"TrackBall"` |
| `selectionMode` | `str` | Selection mode: `"Lasso"` or `"Click"` | `"Lasso"` |

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
    backgroundColor="#000000"
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
    showEnergyPlot=True
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

### Save to File

```python
viewer = MolecularViewer(atoms, style="glossy")
viewer.save_html("molecule_viewer.html")
```
