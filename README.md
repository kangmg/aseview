
<img width="1241" height="598" alt="img" src="https://github.com/user-attachments/assets/701726f6-a23a-48c5-9e1e-34bfb06b8d3a" />

# aseview

> [!WARNING]
> This version (`aseview`) is an alpha release, and detailed instructions will be revised in future updates.
> A new version (v2) is currently under development. It is [pre-alpha](https://github.com/kangmg/aseview_v2_dev) version.
>
> **Known Issue**: Due to a CDN (Content Delivery Network) issue with Java, `aseview` may not function correctly in some environments (e.g., vs-code cursor etc.), excluding Google Colab. This issue will be addressed in a future update.
> This issue was resolved in v2, but further improvements are under development.

`aseview` is a Python package for visualizing atomic structures.

## Installation

You can install `aseview` directly from the GitHub repository using pip:

```bash
pip install git+https://github.com/kangmg/aseview.git -q
```

Alternatively, you can clone the repository and install it locally:

```bash
git clone https://github.com/kangmg/aseview.git
cd aseview
pip install . -q
```

## Sample Data

`aseview` works with `ase.Atoms` objects. You can create simple molecules using `ase.build.molecule` or load structures from various file formats using `ase.io`.

For demonstration and examples, you can download sample trajectory and conformer files using `wget`:

```bash
wget -q https://github.com/kangmg/wget_repo/raw/refs/heads/main/irc_rearraged.traj -O butene.traj
wget -q https://raw.githubusercontent.com/kangmg/wget_repo/refs/heads/main/crest_conformers.xyz -O conformers.xyz
```

## Usage Examples

`aseview` is designed to be used within a Jupyter Notebook or similar environment, as its `viewer` function returns an `IPython.display.HTML` object for interactive 3D visualization.

### Basic Visualization

To visualize an `ase.Atoms` object:

```python
from ase.build import molecule
from aseview import viewer

atoms = molecule('H2O')

viewer(atoms)
```

### Visualizing a Trajectory

You can also visualize a trajectory (a list of `ase.Atoms` objects or an `ase.io.trajectory.TrajectoryReader`). We'll use the `butene.traj` file downloaded in the Sample Data section.

```python
from ase.io import Trajectory
from aseview import viewer

traj = Trajectory('butene.traj')

viewer(traj)
```

### Overlaying Multiple Structures

The `overlay_viewer` function allows you to visualize multiple `ase.Atoms` objects superimposed on each other. This is particularly useful for comparing different conformations or structures. The `superimpose_atoms` function can be used internally to align the structures, or you can manually select specific structures. We'll use the `conformers.xyz` file downloaded in the Sample Data section.

```python
from aseview import overlay_viewer
from ase.io import read

conformers = read('conformers.xyz', index=':', format='xyz')

len(conformers)

first = conformers[0]
first.info['name'] = 'first'

mid = conformers[int(len(conformers) / 2)]
mid.info['name'] = 'mid'

last = conformers[-1]
last.info['name'] = 'last'


overlay_viewer([first, mid, last], opacity=1, bg_color='#000000')
# overlay_viewer(conformers, option='sa', option_param=[0, 1, 2, 4, 6, 10], bg_color='#000000')
```

### Visualizing Normal Modes

The `normal_mode_viewer` function is used to visualize vibrational modes of a molecule. It requires an `ase.Atoms` object and an `ase.vibrations.Vibrations` object.

Here's an example demonstrating how to calculate and visualize normal modes for a benzene molecule:

```python
from ase.build import molecule
from ase.vibrations import Vibrations
from ase.calculators.emt import EMT
from aseview import normal_mode_viewer

atoms = molecule('C6H6')

atoms.calc = EMT()

vib = Vibrations(atoms, delta=0.01, name='benz_vib', nfree=4)
vib.run()

normal_mode_viewer(atoms=atoms, vib=vib, n_frames=20, initial_mode_index=21, bg_color='#ffffff', show_vectors=True, amplitude=2.0, atom_scale=0.6)
```

### Fragment Selection

The `fragment_selector` function provides an interactive viewer for selecting fragments or subsets of atoms from a molecule, with linked 2D and 3D views.

```python
from aseview import fragment_selector
from ase.io import read

bmim = read('conformers.xyz', index='-1', format='xyz')

fragment_selector(bmim)
```
