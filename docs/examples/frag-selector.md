# Fragment Selector

`FragSelector` provides a split 2D/3D interface for interactively choosing atom subsets.
Atoms selected in either panel are kept in sync; the resulting index lists can be copied to the clipboard and used directly in downstream Python code.

---

## Basic Usage

```python
from ase.io import read
from aseview import FragSelector

atoms = read("molecule.xyz")
viewer = FragSelector(atoms)
viewer.show()
```

---

## Selecting Atoms and Using the Indices

Click atoms in either panel to build up a selection, then copy the index list from the sidebar.

```python
from ase.build import molecule
from aseview import FragSelector

ethanol = molecule("CH3CH2OH")
viewer = FragSelector(ethanol)
viewer.show()

# After interacting with the viewer, copy the indices shown in the sidebar:
# e.g. selected   → [0, 1, 2]
# e.g. unselected → [3, 4, 5, 6, 7, 8]

# Use them to slice the Atoms object:
fragment = ethanol[[0, 1, 2]]
```

---

## Working with Fragments (Sub-structures)

A common workflow: visually identify two fragments, copy their indices, then
process each fragment separately.

```python
from ase.io import read
from aseview import FragSelector

# Complex with two ligands
complex_mol = read("complex.xyz")
viewer = FragSelector(complex_mol)
viewer.show()

# After selecting the metal centre + first ligand in the viewer:
metal_ligand1_indices = [0, 3, 5, 7, 8, 9]    # copied from sidebar
ligand2_indices       = [1, 2, 4, 6, 10, 11]   # copied from sidebar

fragment1 = complex_mol[metal_ligand1_indices]
fragment2 = complex_mol[ligand2_indices]
```

---

## Multi-Fragment Molecules

`FragSelector` automatically places disconnected fragments side-by-side in the 2D panel,
making it easy to select entire sub-structures at once.

```python
from ase.io import read
from aseview import FragSelector

# Ion pair, co-crystal, or cluster — fragments appear as separate islands
viewer = FragSelector(read("ion_pair.xyz"))
viewer.show()
```

Use **Select All** → deselect the unwanted fragment → **Copy** to quickly isolate one component.

---

## Adjusting Bond Detection

If bonds are missing or spurious, adjust `bondThreshold` in the sidebar slider,
or set it programmatically:

```python
# Looser threshold — picks up longer / weaker bonds
viewer = FragSelector(atoms, bondThreshold=1.4)
viewer.show()

# Tighter threshold — only very short bonds
viewer = FragSelector(atoms, bondThreshold=1.0)
viewer.show()
```

Changing the slider in the UI rebuilds both the 2D layout and the 3D bonds live.

---

## Custom Appearance

```python
viewer = FragSelector(
    atoms,
    style="glossy",
    atomSize=0.5,
    backgroundColor="#0d1117",
)
viewer.show()
```

---

## Save as Standalone HTML

The generated file is fully self-contained (all JS inlined) and works offline.

```python
viewer = FragSelector(atoms)
viewer.save_html("selector.html")
```

---

## Typical Downstream Workflows

### Constrained Optimisation

```python
from ase.optimize import BFGS
from aseview import FragSelector

atoms = read("structure.xyz")
viewer = FragSelector(atoms)
viewer.show()

# After selecting the frozen core in the viewer:
frozen_indices = [0, 1, 2, 3, 4]   # from sidebar

from ase.constraints import FixAtoms
atoms.set_constraint(FixAtoms(indices=frozen_indices))
opt = BFGS(atoms)
opt.run(fmax=0.05)
```

### QM/MM Region Assignment

```python
# QM region = selected atoms (shown in sidebar as "Selected")
# MM region = unselected atoms (shown as "Unselected")
qm_indices = [0, 1, 2, 5, 6]     # copied from sidebar
mm_indices  = [3, 4, 7, 8, 9]    # copied from sidebar
```

### Trajectory Analysis on a Sub-structure

```python
from ase.io import read
from aseview import FragSelector
import numpy as np

traj   = read("md.traj", index=":")
viewer = FragSelector(traj[0])   # first frame to pick the fragment
viewer.show()

frag_indices = [4, 5, 6, 7]      # selected in viewer

# Compute RMSD of fragment over trajectory
ref_pos = traj[0].get_positions()[frag_indices]
rmsds   = []
for frame in traj:
    pos  = frame.get_positions()[frag_indices]
    diff = pos - ref_pos
    rmsds.append(np.sqrt((diff**2).sum(axis=1).mean()))
```
