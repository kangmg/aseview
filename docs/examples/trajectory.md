# Trajectory Animation

This example demonstrates trajectory visualization with animation and energy plots.

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

- **Play/Pause**: Start/stop animation
- **Frame slider**: Jump to specific frame
- **Speed control**: Adjust playback speed (ms per frame)

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
- Highlights current frame with a marker
- Can be resized by dragging the resize handle

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
