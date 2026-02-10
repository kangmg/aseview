# Normal Mode Visualization

Visualize molecular vibrations from frequency calculations.

## From ORCA Hessian File

=== "Python"

    ```python
    from ase.io import read
    from aseview import NormalViewer

    atoms = read("molecule.xyz")
    viewer = NormalViewer.from_orca(atoms, "orca.hess")
    viewer.show()
    ```

=== "CLI"

    ```bash
    aseview2 molecule.xyz --hess orca.hess
    ```

## From ASE Vibrations

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

## Interactive Controls

### Mode Selection

- **Mode dropdown**: Select vibrational mode by frequency
- **Frequency display**: Shows wavenumber in cm⁻¹
- **Imaginary flag**: Indicates imaginary frequencies (e.g., `123.45i`)

### Animation

- **Play/Pause**: Toggle vibration animation
- **Amplitude slider**: Adjust displacement magnitude

### Mode Vectors

Show atomic displacement arrows:

```python
viewer = NormalViewer.from_orca(
    atoms,
    "orca.hess",
    showModeVector=True
)
```

Vector colors:

- **Red arrows**: Positive displacement direction (+)
- **Blue arrows**: Negative displacement direction (-)

## Customizing Display

### Larger Amplitude

```python
viewer = NormalViewer.from_orca(
    atoms,
    "orca.hess",
    displacementAmplitude=1.5  # Larger vibration
)
```

### Smoother Animation

```python
viewer = NormalViewer(
    atoms,
    vibrations=vib,
    n_frames=60,        # More frames per cycle
    animationSpeed=20   # Faster playback
)
```

### Start at Specific Mode

```python
viewer = NormalViewer.from_orca(
    atoms,
    "orca.hess",
    initialModeIndex=6  # Start at 7th mode (0-indexed)
)
```

## Frequency Types

| Display | Meaning |
|---------|---------|
| `1234.56` | Real frequency (stable vibration) |
| `123.45i` | Imaginary frequency (unstable mode) |

Imaginary frequencies indicate:

- Transition states (one imaginary)
- Higher-order saddle points (multiple imaginary)
- Incomplete optimization

## ORCA Workflow Example

1. **Optimize geometry in ORCA**:
   ```
   ! B3LYP def2-SVP Opt

   * xyz 0 1
   ...
   *
   ```

2. **Run frequency calculation**:
   ```
   ! B3LYP def2-SVP Freq

   * xyzfile 0 1 optimized.xyz
   ```

3. **Visualize in aseview2**:
   ```bash
   aseview2 optimized.xyz --hess orca.hess
   ```

## Supported Formats

| Program | File | Status |
|---------|------|--------|
| ORCA | `.hess` | ✅ Supported |
| VASP | `OUTCAR` | 🔜 Planned |
| Gaussian | `.fchk` | 🔜 Planned |
| xTB | `hessian` | 🔜 Planned |

## Saving Normal Mode Viewer

```python
viewer = NormalViewer.from_orca(atoms, "orca.hess")
viewer.save_html("vibrations.html")
```

The saved HTML file includes all vibrational modes and can be viewed without Python.
