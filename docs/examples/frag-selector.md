# Fragment Selector

`FragSelector` provides a split 2D/3D interface for interactively choosing atom subsets.
Atoms selected in either panel are kept in sync; the resulting index lists can be copied to the clipboard and used directly in downstream Python code.

---

## Live Demo

<iframe src="../../assets/viewers/frag_selector_basic.html" width="100%" height="550" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

Click any atom in either panel to toggle its selection. The **Selected** and **Unselected** index lists in the sidebar update live. Hit **Copy** to paste indices into your next cell.

---

## Basic Usage

=== "Python"

    ```python
    from ase.build import molecule
    from aseview import FragSelector

    atoms = molecule('CH3OH')  # methanol
    viewer = FragSelector(atoms)
    viewer.show()
    ```

=== "CLI"

    ```bash
    aseview molecule.xyz -v frag
    ```

---

## Selecting a Fragment and Slicing the Structure

```python
from ase.build import molecule
from aseview import FragSelector

ethanol = molecule('CH3CH2OH')  # 9 atoms: C, C, O, 6×H
viewer = FragSelector(ethanol)
viewer.show()
```

After selecting the hydroxyl group (O + H) in the viewer the sidebar shows, for example:

```
Selected   2  →  [2, 8]
Unselected 7  →  [0, 1, 3, 4, 5, 6, 7]
```

Use those indices directly in your code:

```python
oh_indices    = [2, 8]                   # copied from "Selected"
ethyl_indices = [0, 1, 3, 4, 5, 6, 7]  # copied from "Unselected"

oh_group = ethanol[oh_indices]
ethyl    = ethanol[ethyl_indices]

print(oh_group.get_chemical_formula())   # HO
print(ethyl.get_chemical_formula())      # C2H5
```

---

## Multi-Fragment Molecules

Disconnected components automatically appear as separate islands in the 2D panel,
making it easy to pick an entire sub-structure with a few clicks.

<iframe src="../../assets/viewers/frag_selector_multifrag.html" width="100%" height="550" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

```python
from ase.build import molecule
from aseview import FragSelector

# Build two separated water molecules
water1 = molecule('H2O')
water2 = molecule('H2O')
water2.translate([5, 0, 0])   # move second molecule 5 Å away

system = water1 + water2      # 6 atoms total
viewer = FragSelector(system)
viewer.show()
```

Use **Select All** then click the atoms you *don't* want, or click each molecule's atoms individually.

```python
# After selecting only the first water molecule in the viewer:
mol1_indices = [0, 1, 2]   # from sidebar "Selected"
mol2_indices = [3, 4, 5]   # from sidebar "Unselected"

first_water  = system[mol1_indices]
second_water = system[mol2_indices]
```

---

## Selection Modes

The 2D panel has four interaction modes switchable via the toolbar (top-right of the panel) or keyboard shortcuts:

| Mode | Shortcut | Description |
|------|----------|-------------|
| **Navigate** | <kbd>N</kbd> | Drag to pan, scroll to zoom, double-click to fit |
| **Click** | <kbd>C</kbd> | Click individual atoms to toggle selection |
| **Rect** | <kbd>R</kbd> | Drag a rectangle — all atoms inside are selected |
| **Lasso** | <kbd>L</kbd> | Draw a freeform region — all enclosed atoms are selected |

Hold **Shift** while releasing the mouse to *add* the new atoms to the existing selection.

### Large Conjugated System — Bis-iodo Heterocycles

The viewer below shows two disconnected iodo-substituted benzimidazole ring systems (50 atoms total).
Try switching to **Rect** or **Lasso** mode to select one ring system at a time.

<iframe src="../../assets/viewers/frag_selector_bis_iodo.html" width="100%" height="550" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

```python
from ase import Atoms
from aseview import FragSelector

symbols = [
    'C','C','C','C','C','C','N','C','C','C','C','C','C',
    'H','H','H','H','H','H','H','H','H','I','H','H',
    'C','C','C','C','C','C','N','C','C','C','C','C','C',
    'H','H','H','H','H','H','H','H','H','I','H','H',
]
positions = [
    [ 0.0947,-1.8058, 0.6052], [-0.1026,-0.4427, 0.4521],
    [-0.3876, 0.3096, 1.5873], [-0.4896,-0.3395, 2.8462],
    [-0.3517,-1.7140, 2.9702], [-0.0383,-2.4448, 1.8379],
    [-0.6878, 0.6032, 3.8288], [-0.7234, 1.8569, 3.2472],
    [-0.5405, 1.7201, 1.8456], [-0.4876, 2.8627, 1.0532],
    [-0.6277, 4.1064, 1.6489], [-0.8218, 4.2216, 3.0286],
    [-0.8675, 3.1015, 3.8426], [-0.7044, 0.3061, 5.2390],
    [-0.0290, 0.0295,-0.5177], [-0.4422,-2.2053, 3.9290],
    [ 0.1302,-3.5082, 1.9182], [-0.3359, 2.7830,-0.0158],
    [-0.9305, 5.2037, 3.4696], [-0.9763, 3.2031, 4.9124],
    [ 0.0039,-0.4918, 5.4470], [-0.3875, 1.1795, 5.8004],
    [ 0.7151,-2.9701,-1.0540], [-1.6959,-0.0031, 5.5734],
    [-0.5902, 4.9991, 1.0417],
    [ 2.5237, 2.5761, 5.5852], [ 2.6805, 1.2022, 5.5217],
    [ 2.8972, 0.6006, 4.2861], [ 2.9543, 1.4108, 3.1229],
    [ 2.8109, 2.7881, 3.1876], [ 2.5916, 3.3586, 4.4296],
    [ 3.1388, 0.6022, 2.0188], [ 3.2099,-0.7112, 2.4327],
    [ 3.0517,-0.7623, 3.8435], [ 3.0014,-1.9949, 4.4854],
    [ 3.1289,-3.1417, 3.7182], [ 3.3069,-3.0879, 2.3329],
    [ 3.3639,-1.8689, 1.6816], [ 3.1678, 1.0657, 0.6548],
    [ 2.6368, 0.5998, 6.4199], [ 2.8409, 3.3943, 2.2941],
    [ 2.4565, 4.4288, 4.5049], [ 2.8469,-2.0527, 5.5535],
    [ 3.3795,-3.9994, 1.7582], [ 3.4966,-1.8333, 0.6099],
    [ 4.1569, 1.4375, 0.3810], [ 2.8930, 0.2545,-0.0129],
    [ 2.4368, 1.8602, 0.5232], [ 3.0171,-5.0395, 4.6479],
    [ 2.3587, 3.0522, 6.5411],
]

atoms = Atoms(symbols=symbols, positions=positions)
viewer = FragSelector(atoms, style='cartoon')
viewer.show()
```

Select one ring system using **Lasso** mode, then copy the indices to use them in further computations.

---

## Constrained Geometry Optimisation

Visually pick the atoms you want to freeze, then pass the indices to `FixAtoms`.

```python
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from aseview import FragSelector

ethane = molecule('C2H6')   # 2 carbons + 6 hydrogens
viewer = FragSelector(ethane)
viewer.show()
```

Select one methyl group in the viewer; the sidebar shows something like:

```
Selected   4  →  [0, 2, 3, 4]
```

Freeze those atoms and optimise:

```python
frozen_indices = [0, 2, 3, 4]   # from sidebar "Selected"

ethane.calc = EMT()
ethane.set_constraint(FixAtoms(indices=frozen_indices))
opt = BFGS(ethane, logfile=None)
opt.run(fmax=0.05)

print("Optimisation done. Frozen:", frozen_indices)
```

---

## QM/MM Region Assignment

`FragSelector` maps naturally onto QM/MM workflows: the **Selected** list becomes the QM region and the **Unselected** list becomes the MM region.

```python
from ase.build import molecule
from aseview import FragSelector

# Methanol: C(0) O(1) H×4 = 6 atoms
atoms = molecule('CH3OH')
viewer = FragSelector(atoms)
viewer.show()
```

Select the reactive OH group; the sidebar gives:

```
Selected   2  →  [1, 5]       ← QM region (O + H)
Unselected 4  →  [0, 2, 3, 4] ← MM region (C + 3H)
```

```python
qm_indices = [1, 5]        # from sidebar "Selected"
mm_indices  = [0, 2, 3, 4] # from sidebar "Unselected"

qm_atoms = atoms[qm_indices]
mm_atoms = atoms[mm_indices]
```

---

## Adjusting Bond Detection

Increase or decrease `bondThreshold` if bonds are missing or spurious.
The slider in the **Display** card applies the change instantly — both the 2D layout and 3D bonds rebuild live.

=== "Python"

    ```python
    from ase.io import read
    from aseview import FragSelector

    # Loose threshold: include weak / elongated bonds
    viewer = FragSelector(read("structure.xyz"), bondThreshold=1.4)
    viewer.show()

    # Tight threshold: only very short covalent bonds
    viewer = FragSelector(read("structure.xyz"), bondThreshold=1.0)
    viewer.show()
    ```

=== "CLI"

    ```bash
    # Loose threshold
    aseview structure.xyz -v frag --bond-threshold 1.4

    # Tight threshold
    aseview structure.xyz -v frag --bond-threshold 1.0
    # short form: -bt
    aseview structure.xyz -v frag -bt 1.0
    ```

---

## Custom Appearance

```python
from ase.build import molecule
from aseview import FragSelector

atoms = molecule('C2H6')
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

The generated file bundles all JavaScript inline and works in any browser without Python.

=== "Python"

    ```python
    from ase.build import molecule
    from aseview import FragSelector

    viewer = FragSelector(molecule('CH3OH'))
    viewer.save_html("selector.html")
    ```

=== "CLI"

    ```bash
    aseview molecule.xyz -v frag -o selector.html
    ```

---

## Trajectory: Pick a Sub-structure, Then Analyse It

Use the first frame to choose a fragment visually, then apply the indices across the whole trajectory.

```python
from ase.io import read
from aseview import FragSelector
import numpy as np

traj = read("md.traj", index=":")

# Open first frame to pick the fragment of interest
FragSelector(traj[0]).show()
```

After selecting the fragment in the viewer:

```python
frag_indices = [4, 5, 6, 7]   # from sidebar "Selected"

# Compute per-frame RMSD of the fragment relative to frame 0
ref = traj[0].get_positions()[frag_indices]
rmsds = []
for frame in traj:
    pos = frame.get_positions()[frag_indices]
    rmsds.append(np.sqrt(((pos - ref)**2).sum(axis=1).mean()))
```
