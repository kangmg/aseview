# OverlayViewer

Overlay viewer for comparing multiple molecular structures simultaneously.

## Constructor

```python
OverlayViewer(data, index_list=None, **kwargs)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `data` | `Atoms`, `List[Atoms]`, `Dict`, `str` | One or more molecular structures | Required |
| `index_list` | `None`, `list[int]`, `list[list[int]]` | Atom selection — see table below | `None` |

#### index_list Modes

The mode is inferred automatically from the type of `index_list`:

| Value | Mode | Behaviour |
|-------|------|-----------|
| `None` | **all** | All atoms from every structure are used (default) |
| `[0, 1, 2]` | **same** | Same atom indices applied to every structure |
| `[[0, 1], [2, 3], [0, 2]]` | **indices** | Different atom indices per structure (length must equal number of structures) |

### Keyword Arguments (Settings)

#### Geometry Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atomSize` | `float` | Atom sphere radius scale | `0.4` |
| `bondThickness` | `float` | Bond cylinder radius | `0.1` |
| `bondThreshold` | `float` | Bond detection threshold | `1.2` |

#### Style Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `style` | `str` | Visual style name | `"cartoon"` |
| `backgroundColor` | `str` | Background color (hex) | `"#1f2937"` |

#### Color Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `colorBy` | `str` | Coloring mode | `"Atom"` |
| `colormap` | `str` | Colormap name (when `colorBy="Colormap"`) | `"viridis"` |

##### colorBy Options

| Value | Description |
|-------|-------------|
| `"Atom"` | Color by element (CPK) |
| `"Molecule"` | Each molecule has distinct color |
| `"Colormap"` | Gradient colormap based on molecule index |

##### Available Colormaps

| Name | Description |
|------|-------------|
| `"viridis"` | Perceptually uniform (blue → green → yellow) |
| `"plasma"` | Perceptually uniform (purple → orange → yellow) |
| `"coolwarm"` | Diverging (blue → white → red) |
| `"jet"` | Rainbow (blue → cyan → yellow → red) |
| `"rainbow"` | Full spectrum |
| `"grayscale"` | Black to white |

#### Display Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `showCell` | `bool` | Show unit cell | `True` |
| `showBond` | `bool` | Show bonds | `True` |
| `showShadow` | `bool` | Enable shadows | `False` |
| `centerMolecules` | `bool` | Center each molecule at origin before overlay | `False` |

#### Alignment Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `alignMolecules` | `bool` | Align molecules using Kabsch + Hungarian algorithm | `False` |

!!! tip "Molecular Alignment"
    When `alignMolecules=True`:

    - First molecule is used as reference
    - Other molecules are aligned using RMSD minimization
    - **Kabsch algorithm**: Finds optimal rotation matrix
    - **Hungarian reordering**: Matches atoms by element type and position
    - Automatically centers molecules (no need for `centerMolecules`)

#### View Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `viewMode` | `str` | `"Perspective"` or `"Orthographic"` | `"Perspective"` |
| `rotationMode` | `str` | `"TrackBall"` or `"Orbit"` | `"TrackBall"` |

## Methods

### show()

Display the viewer in a Jupyter notebook.

```python
viewer.show(width='100%', height=600)
```

### get_html()

Get the HTML content as a string.

```python
html = viewer.get_html()
```

### save_html()

Save the viewer as an HTML file.

```python
viewer.save_html(filename)
```

## UI Controls (Interactive)

The overlay viewer provides interactive controls:

### Molecule Panel

- **Visibility toggle**: Show/hide individual molecules
- **Opacity slider**: Adjust transparency (0-100%)
- **Color picker**: Change molecule color

### Animation Controls

- **Play/Pause**: Animate through frames
- **Frame slider**: Manual frame selection
- **Speed control**: Adjust playback speed

## Examples

### Compare Two Structures

```python
from ase.io import read
from aseview import OverlayViewer

reactant = read("reactant.xyz")
product = read("product.xyz")

viewer = OverlayViewer([reactant, product])
viewer.show()
```

### Atom Selection with index_list

```python
# all (default) – show every atom
viewer = OverlayViewer([mol1, mol2, mol3])

# same – show atoms 0, 1, 2 from every structure
viewer = OverlayViewer([mol1, mol2, mol3], index_list=[0, 1, 2])

# indices – different atoms per structure
viewer = OverlayViewer(
    [mol1, mol2, mol3],
    index_list=[[0, 1, 2], [0, 1, 3], [1, 2, 3]]
)
viewer.show()
```

### Trajectory with Colormap

```python
trajectory = read("optimization.xyz", index=":")

viewer = OverlayViewer(
    trajectory,
    colorBy="Colormap",
    colormap="viridis"
)
viewer.show()
```

### Centered Molecules

```python
# Center each molecule at origin for better comparison
viewer = OverlayViewer(
    [mol1, mol2, mol3],
    centerMolecules=True,
    colorBy="Molecule"
)
viewer.show()
```

### Custom Molecule Names

Set custom names for molecules using `atoms.info['name']`:

```python
from ase.build import molecule
from aseview import OverlayViewer

# Create molecules and set names
reactant = molecule("C2H6")
reactant.info['name'] = "Ethane (reactant)"

ts = molecule("C2H6")
ts.info['name'] = "Transition State"

product = molecule("C2H4")
product.info['name'] = "Ethene (product)"

viewer = OverlayViewer(
    [reactant, ts, product],
    colorBy="Molecule",
    centerMolecules=True
)
viewer.show()
```

!!! tip "Naming Convention"
    Without `atoms.info['name']`, molecules are labeled as "Molecule 1", "Molecule 2", etc.
    Custom names appear in the molecule control panel for easier identification.

### Aligned Molecules (RMSD Minimization)

```python
# Align molecules using Kabsch rotation + Hungarian reordering
# First molecule is the reference
viewer = OverlayViewer(
    [reference, conformer1, conformer2],
    alignMolecules=True,
    colorBy="Molecule"
)
viewer.show()
```

This uses:

- **Kabsch algorithm**: Optimal rotation via SVD
- **Hungarian algorithm**: Atom correspondence matching by element type

### Custom Styling

```python
viewer = OverlayViewer(
    structures,
    style="glossy",
    colorBy="Colormap",
    colormap="plasma",
    atomSize=0.5
)
viewer.show()
```

### Save Comparison

```python
viewer = OverlayViewer([before, after], colorBy="Molecule")
viewer.save_html("comparison.html")
```
