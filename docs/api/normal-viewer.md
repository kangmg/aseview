# NormalViewer

Vibrational normal mode viewer for molecular frequency analysis.

## Constructor

```python
NormalViewer(atoms, vibrations=None, mode_vectors=None, frequencies=None, n_frames=30, **kwargs)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atoms` | `Atoms`, `Dict` | Equilibrium molecular structure | Required |
| `vibrations` | `Vibrations`, `VibrationsData` | ASE vibrations object | `None` |
| `mode_vectors` | `List` | Mode displacement vectors | `None` |
| `frequencies` | `List` | Frequencies in cm⁻¹ | `None` |
| `n_frames` | `int` | Animation frames per cycle | `30` |

!!! note
    Provide either `vibrations` OR both `mode_vectors` and `frequencies`.

### Keyword Arguments (Settings)

#### Geometry Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atomSize` | `float` | Atom sphere radius scale | `0.4` |
| `bondThickness` | `float` | Bond cylinder radius | `0.1` |
| `bondThreshold` | `float` | Bond detection threshold | `1.0` |

#### Style Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `style` | `str` | Visual style name | `"cartoon"` |
| `backgroundColor` | `str` | Background color (hex) | `"#1f2937"` |

#### Display Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `showCell` | `bool` | Show unit cell | `True` |
| `showBond` | `bool` | Show bonds | `True` |
| `showShadow` | `bool` | Enable shadows | `False` |
| `showModeVector` | `bool` | Show displacement arrows | `False` |

#### Animation Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `animationSpeed` | `int` | Animation speed (ms) | `30` |
| `displacementAmplitude` | `float` | Vibration amplitude scale | `0.75` |
| `initialModeIndex` | `int` | Starting mode index | `0` |
| `nFrames` | `int` | Frames per vibration cycle | `30` |

## Class Methods

### from_orca()

Create NormalViewer from ORCA Hessian file.

```python
NormalViewer.from_orca(atoms, hess_file, skip_imaginary=False, **kwargs)
```

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atoms` | `Atoms` | Equilibrium structure | Required |
| `hess_file` | `str` | Path to ORCA `.hess` file | Required |
| `skip_imaginary` | `bool` | Filter out modes with freq < 10 cm⁻¹ | `False` |

## Instance Methods

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

### Normal Modes Panel

- **Mode selector**: Choose vibration mode by frequency
- **Frequency display**: Shows current mode frequency (cm⁻¹)
- **Imaginary indicator**: Highlights imaginary frequencies

### Animation Controls

- **Play/Pause**: Toggle vibration animation
- **Amplitude slider**: Adjust displacement magnitude
- **Mode vectors toggle**: Show/hide displacement arrows

### Mode Vector Display

When enabled, arrows show atomic displacements:

- **Red arrows**: Positive direction (+)
- **Blue arrows**: Negative direction (-)

## Examples

### From ASE Vibrations

```python
from ase.vibrations import Vibrations
from aseview import NormalViewer

# Run frequency calculation
vib = Vibrations(atoms)
vib.run()

# Visualize
viewer = NormalViewer(atoms, vibrations=vib)
viewer.show()
```

### From ORCA Hessian

```python
from ase.io import read
from aseview import NormalViewer

atoms = read("molecule.xyz")
viewer = NormalViewer.from_orca(atoms, "orca.hess")
viewer.show()
```

### With Mode Vectors Visible

```python
viewer = NormalViewer.from_orca(
    atoms,
    "orca.hess",
    showModeVector=True,
    displacementAmplitude=1.0
)
viewer.show()
```

### Custom Animation

```python
viewer = NormalViewer(
    atoms,
    vibrations=vib,
    n_frames=60,           # Smoother animation
    animationSpeed=20,     # Faster playback
    initialModeIndex=6     # Start at 7th mode
)
viewer.show()
```

### Save Normal Mode Viewer

```python
viewer = NormalViewer.from_orca(atoms, "orca.hess")
viewer.save_html("vibrations.html")
```

## Supported Hessian Formats

| Program | File | Status |
|---------|------|--------|
| ORCA | `.hess` | Supported |
| VASP | `OUTCAR` | Planned |
| Gaussian | `.fchk` | Planned |
| xTB | `hessian` | Planned |

## Frequency Display

Frequencies are displayed with special formatting:

| Type | Display | Meaning |
|------|---------|---------|
| Real | `1234.56` | Normal vibration |
| Imaginary | `123.45i` | Transition state / unstable mode |
