# Normal Mode Visualization

Visualize molecular vibrations from frequency calculations.

## Live Demo

<iframe src="../../assets/viewers/normal_modes.html" width="100%" height="600" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## From ORCA Hessian

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

---

## From ASE Vibrations

```python
from ase.vibrations import Vibrations
from aseview import NormalViewer

vib = Vibrations(atoms)
vib.run()

viewer = NormalViewer(atoms, vibrations=vib)
viewer.show()
```

---

## Options

```python
viewer = NormalViewer.from_orca(
    atoms,
    "orca.hess",
    showModeVector=True,        # Show displacement arrows
    displacementAmplitude=1.5,  # Larger vibration
    initialModeIndex=6          # Start at 7th mode
)
```

---

## Interactive Controls

| Control | Description |
|---------|-------------|
| **Mode dropdown** | Select vibrational mode |
| **Play/Pause** | Toggle animation |
| **Amplitude slider** | Adjust displacement |
| **Mode Vectors** | Show/hide displacement arrows |

---

## Frequency Display

| Format | Meaning |
|--------|---------|
| `1234.56` | Real frequency (stable) |
| `123.45i` | Imaginary frequency (unstable) |

Imaginary frequencies indicate transition states or incomplete optimization.

---

## Supported Formats

| Program | File | Status |
|---------|------|--------|
| ORCA | `.hess` | Supported |
| VASP | `OUTCAR` | Planned |
| Gaussian | `.fchk` | Planned |

---

## Saving

```python
viewer = NormalViewer.from_orca(atoms, "orca.hess")
viewer.save_html("vibrations.html")
```
