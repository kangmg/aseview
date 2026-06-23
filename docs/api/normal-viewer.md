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
| `theme` | `str` | Visual theme (`"dark"`, `"spring"`, `"glass"`, `"darkgreen"`, `"simple"`, …). `None` uses the global default. | `None` |

!!! note
    Provide either `vibrations` OR both `mode_vectors` and `frequencies`.

### Keyword Arguments (Settings)

#### Geometry Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atomSize` | `float` | Atom sphere radius scale | `0.4` |
| `bondThickness` | `float` | Bond cylinder radius | `0.09` |
| `bondThreshold` | `float` | Bond detection threshold | `1.2` |

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
| `hideHydrogens` | `bool` | Hide hydrogen atoms and their bonds without removing them from the data | `False` |
| `showBlur` | `bool` | Apply blur to the rendered molecule canvas | `False` |
| `blurStrength` | `float` | Blur radius in pixels when `showBlur` is enabled | `1.5` |
| `showShadow` | `bool` | Enable shadows | `False` |
| `showModeVector` | `bool` | Show displacement arrows | `False` |

#### Animation Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `animationSpeed` | `int` | Animation speed (ms) | `30` |
| `displacementAmplitude` | `float` | Vibration amplitude scale | `0.75` |
| `initialModeIndex` | `int` | Starting mode index | `0` |
| `nFrames` | `int` | Frames per vibration cycle | `30` |

#### View Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `viewMode` | `str` | Camera projection: `"Orthographic"` or `"Perspective"` | `"Orthographic"` |
| `viewPreset` | `str`, `None` | Initial named camera direction (`"top-c"`, `"side-a"`, `"front"`, etc.) | `None` |
| `viewDirection` | `list[float]`, `None` | Explicit target-to-camera direction vector | `None` |
| `viewEuler` | `list[float]`, `None` | `[rx, ry, rz]` degrees in XYZ order, applied to `[0, 0, 1]` | `None` |
| `viewUp` | `list[float]`, `None` | Optional camera up-vector hint | `None` |
| `viewFit` | `float` | Camera fit multiplier | `1.0` |

`viewDirection` and preset directions are target-to-camera vectors. Presets
`top`, `bottom`, `front`, `back`, `left`, and `right` use Cartesian axes.
Presets `top-c`, `bottom-c`, `side-a`, and `side-b` use valid unit-cell
vectors. Aliases `c` and `top` prefer the cell `c` axis, while `a` and `b`
prefer the cell `a` and `b` axes. Missing, non-periodic, zero-length, or
degenerate cells fall back to Cartesian directions.

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

### from_vasp()

Create NormalViewer from VASP OUTCAR file (IBRION=5 or 6).

```python
NormalViewer.from_vasp(atoms, outcar_file, skip_imaginary=False, **kwargs)
```

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atoms` | `Atoms` | Equilibrium structure | Required |
| `outcar_file` | `str` | Path to VASP OUTCAR file | Required |
| `skip_imaginary` | `bool` | Filter out modes with freq < 10 cm⁻¹ | `False` |

!!! note
    The OUTCAR must contain the "Eigenvectors and eigenvalues of the dynamical matrix" section produced by `IBRION=5` or `IBRION=6`. Imaginary (soft) modes are reported with negative frequencies.

### from_file()

Auto-detect format and create NormalViewer from any supported Hessian file.

```python
NormalViewer.from_file(atoms, filepath, skip_imaginary=False, **kwargs)
```

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `atoms` | `Atoms` | Equilibrium structure | Required |
| `filepath` | `str` | Path to ORCA `.hess` or VASP `OUTCAR` | Required |
| `skip_imaginary` | `bool` | Filter out modes with freq < 10 cm⁻¹ | `False` |

Format detection order: filename (`OUTCAR` → VASP), then content sniffing.

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

!!! note "Image export"
    `NormalViewer` can save PNG/GIF images from Python when the optional export
    dependency is installed:

    ```bash
    pip install "aseview[export]"
    python -m playwright install chromium
    ```

    ```python
    viewer.save_png("normal-mode.png", scale=2)
    viewer.save_gif("normal-mode.gif", frames=30, quality="high")
    ```

    PNG quality is controlled by `scale`, `width`, and `height`; GIF
    `quality="low" | "medium" | "high"` maps to encoder sampling quality.
    The CLI does not provide `--save-png` or `--save-gif`.

    In the browser JavaScript module, `ASEView.NormalModeViewer` exposes
    `setView(viewSpec)`, `resetView()`, `savePNG(options)`, and
    `saveGIF(options)`. Normal GIF export animates the current normal mode.
    Export options are `filename`, `download`, `returnDataUrl`, `scale`,
    `width`, `height`, `transparent`, `backgroundColor`, plus GIF-only
    `frames`, `delay`, and `sampleInterval`. Export promises resolve with
    `{ ok: true, type: "png" | "gif", filename, dataUrl?, width?, height? }` or
    reject with an `Error` carrying `code`, `message`, `type`, and when
    available `requestId`.

## UI Controls (Interactive)

### Normal Modes Panel

- **Mode selector**: Choose vibration mode by frequency
- **Frequency display**: Shows current mode frequency (cm⁻¹)
- **Imaginary indicator**: Highlights imaginary frequencies

### Animation Controls

- **Play/Pause**: Toggle vibration animation
- **Amplitude slider**: Adjust displacement magnitude
- **Mode vectors toggle**: Show/hide displacement arrows
- **Copy format chip**: Cycle structure copy format: `xyz`, `extxyz`, `cif`, `POSCAR`
- **Copy frame / Copy all**: Copy the current displaced structure or all generated vibration frames

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

### From VASP OUTCAR

```python
from ase.io import read
from aseview import NormalViewer

atoms = read("POSCAR")
viewer = NormalViewer.from_vasp(atoms, "OUTCAR")
viewer.show()
```

### Auto-detect format

```python
from ase.io import read
from aseview import NormalViewer

atoms = read("molecule.xyz")
viewer = NormalViewer.from_file(atoms, "orca.hess")   # ORCA
viewer = NormalViewer.from_file(atoms, "OUTCAR")       # VASP
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
viewer = NormalViewer.from_file(atoms, "orca.hess")
viewer.save_html("vibrations.html")
```

## Supported Hessian Formats

| Program | File | Status |
|---------|------|--------|
| ORCA | `.hess` | Supported |
| VASP | `OUTCAR` (IBRION=5/6) | Supported |
| Gaussian | `.fchk` | Planned |
| xTB | `hessian` | Planned |

## Frequency Display

Frequencies are displayed with special formatting:

| Type | Display | Meaning |
|------|---------|---------|
| Real | `1234.56` | Normal vibration |
| Imaginary | `123.45i` | Transition state / unstable mode |
