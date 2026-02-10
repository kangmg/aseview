# Trajectory Animation

This example demonstrates trajectory visualization with animation and energy plots.

## Live Demo: Trajectory with Energy Plot

<iframe src="../../assets/viewers/trajectory.html" width="100%" height="600" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Loading a Trajectory

=== "Python"

    ```python
    from ase.io import read
    from aseview import MolecularViewer

    # Load all frames from trajectory file
    trajectory = read("md_trajectory.xyz", index=":")

    viewer = MolecularViewer(trajectory)
    viewer.show()
    ```

=== "CLI"

    ```bash
    # All frames
    aseview2 md_trajectory.xyz

    # Subset of frames
    aseview2 md_trajectory.xyz -i 0:100

    # Every 10th frame
    aseview2 md_trajectory.xyz -i ::10
    ```

## Animation Controls

The viewer provides playback controls:

| Control | Description |
|---------|-------------|
| **Play/Pause** | Start/stop animation |
| **Frame slider** | Jump to specific frame |
| **Speed control** | Adjust playback speed (ms per frame) |
| **First/Last** | Jump to first or last frame |

## Energy Plot

If energy data is available in the trajectory, you can display an energy plot synchronized with the animation:

```python
viewer = MolecularViewer(
    trajectory,
    showEnergyPlot=True
)
viewer.show()
```

The energy plot:

- Shows energy vs. frame number
- Highlights current frame with a red marker
- Resizable by dragging the handle at the bottom

### Energy Data Sources

Energy is automatically extracted from:

1. **ASE Calculator**: `atoms.get_potential_energy()`
2. **Info dict**: `atoms.info['energy']`, `atoms.info['Energy']`, etc.

## Force Vectors

Display force vectors on atoms (if force data is available):

```python
viewer = MolecularViewer(
    trajectory,
    showForces=True,
    forceScale=2.0  # Adjust arrow length
)
viewer.show()
```

## Frame Selection Examples

| Index | Result |
|-------|--------|
| `-i 0` | First frame only |
| `-i -1` | Last frame only |
| `-i :` | All frames |
| `-i 0:10` | Frames 0-9 |
| `-i ::2` | Every 2nd frame |
| `-i 10:50:5` | Frames 10,15,20,...,45 |
| `-i -20:` | Last 20 frames |

## Optimization Trajectory

Visualizing geometry optimization:

```python
# Load optimization trajectory
opt_traj = read("optimization.traj", index=":")

viewer = MolecularViewer(
    opt_traj,
    showEnergyPlot=True,  # See energy decrease
    animationSpeed=100     # Slower playback
)
viewer.show()
```

## Saving Trajectory Viewer

```python
viewer = MolecularViewer(trajectory, showEnergyPlot=True)
viewer.save_html("trajectory_viewer.html")
```

The saved HTML includes all frames and can be animated in the browser.
