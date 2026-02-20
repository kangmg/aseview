# Trajectory Animation

Visualize molecular dynamics trajectories with animation and energy plots.

## Live Demo

<iframe src="../../assets/viewers/trajectory.html" width="100%" height="600" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Basic Usage

=== "Python"

    ```python
    from ase.io import read
    from aseview import MolecularViewer

    trajectory = read("md_trajectory.xyz", index=":")
    viewer = MolecularViewer(trajectory)
    viewer.show()
    ```

=== "CLI"

    ```bash
    aseview md_trajectory.xyz
    ```

---

## Frame Selection

| Index | Result |
|-------|--------|
| `-i :` | All frames |
| `-i 0:10` | Frames 0-9 |
| `-i ::2` | Every 2nd frame |
| `-i -1` | Last frame only |
| `-i -20:` | Last 20 frames |

---

## Energy Plot

Display energy synchronized with animation:

```python
viewer = MolecularViewer(trajectory, showEnergyPlot=True)
viewer.show()
```

Energy is extracted from `atoms.get_potential_energy()` or `atoms.info['energy']`.

---

## Force Vectors

Display force vectors on atoms:

```python
viewer = MolecularViewer(trajectory, showForces=True, forceScale=2.0)
viewer.show()
```

---

## Animation Controls

| Control | Description |
|---------|-------------|
| **Play/Pause** | Toggle animation |
| **Frame slider** | Jump to frame |
| **Copy frame** | Copy current frame as XYZ |
| **Copy all** | Copy entire trajectory |

---

## Saving

```python
viewer = MolecularViewer(trajectory, showEnergyPlot=True)
viewer.save_html("trajectory.html")
```
