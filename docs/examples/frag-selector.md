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

Both the **2D** and **3D** panels share a single mode toolbar (centred at the top). Switch modes via the toolbar buttons or keyboard shortcuts:

| Mode | Shortcut | 2D behaviour | 3D behaviour |
|------|----------|--------------|--------------|
| **Navigate** | <kbd>N</kbd> | Drag to pan, scroll to zoom, double-click to fit | Drag to rotate, right-drag to pan, scroll to zoom |
| **Click** | <kbd>C</kbd> | Click atom to toggle selection | Click atom to toggle selection, drag still rotates |
| **Rect** | <kbd>R</kbd> | Drag a rectangle — all enclosed atoms selected | Drag a rectangle on-screen — all atoms whose 3D position projects inside are selected |
| **Lasso** | <kbd>L</kbd> | Draw a freeform region — all enclosed atoms selected | Draw a freeform region on-screen — all atoms whose 3D position projects inside are selected |

Hold **Shift** while releasing the mouse to *add* the new atoms to the existing selection.

### Large Conjugated System — Bis-iodo Heterocycles

The viewer below shows two disconnected iodo-substituted benzimidazole ring systems (50 atoms total).
Try switching to **Rect** or **Lasso** mode and drawing a box around one ring system in either the 2D or 3D panel.

<iframe src="../../assets/viewers/frag_selector_bis_iodo.html" width="100%" height="550" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

```python
from ase.io import read
from aseview import FragSelector

atoms = read("dimer.xyz")   # two iodo-benzimidazole units, 50 atoms
viewer = FragSelector(atoms, style='cartoon')
viewer.show()
```

Select one ring system using **Lasso** mode in the 2D or 3D panel, then copy the indices to use them in further computations.

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
