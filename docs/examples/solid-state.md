# Solid-State Structures

Visualize crystals, surfaces, and periodic systems with unit cell display.

## Live Demo: Silicon Crystal

Silicon in diamond structure (2x2x2 supercell):

<iframe src="../../assets/viewers/silicon_crystal.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Live Demo: NaCl Crystal

Sodium chloride in rocksalt structure:

<iframe src="../../assets/viewers/nacl_crystal.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Live Demo: Au(111) Surface

Gold surface slab (4 layers, 3x3 cell):

<iframe src="../../assets/viewers/gold_surface.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Live Demo: Graphene

Graphene nanoribbon (zigzag edge):

<iframe src="../../assets/viewers/graphene.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Live Demo: Carbon Nanotube

(6,0) Carbon nanotube:

<iframe src="../../assets/viewers/carbon_nanotube.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Live Demo: FCC Copper

Copper FCC crystal (3x3x3 supercell):

<iframe src="../../assets/viewers/copper_fcc.html" width="100%" height="500" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

## Building Crystals

### Bulk Structures

=== "Python"

    ```python
    from ase.build import bulk
    from aseview import MolecularViewer

    # Diamond structure (Si, C, Ge)
    si = bulk('Si', 'diamond', a=5.43, cubic=True)
    si = si * (2, 2, 2)  # 2x2x2 supercell

    viewer = MolecularViewer(si, showCell=True)
    viewer.show()
    ```

=== "CLI"

    ```bash
    aseview2 POSCAR --style metallic
    aseview2 structure.cif
    ```

### Common Crystal Structures

| Structure | ASE Function | Example |
|-----------|-------------|---------|
| FCC | `bulk('Cu', 'fcc', a=3.61)` | Cu, Ag, Au, Al, Ni |
| BCC | `bulk('Fe', 'bcc', a=2.87)` | Fe, W, Cr, Mo |
| Diamond | `bulk('Si', 'diamond', a=5.43)` | Si, Ge, C |
| Rocksalt | `bulk('NaCl', 'rocksalt', a=5.64)` | NaCl, MgO, LiF |
| Zincblende | `bulk('GaAs', 'zincblende', a=5.65)` | GaAs, ZnS |
| Wurtzite | `bulk('ZnO', 'wurtzite', ...)` | ZnO, GaN |
| HCP | `bulk('Mg', 'hcp', a=3.21, c=5.21)` | Mg, Ti, Zn |

### Supercells

```python
# Create supercell
atoms = bulk('Si', 'diamond', a=5.43)
supercell = atoms * (3, 3, 3)  # 3x3x3 supercell

viewer = MolecularViewer(supercell, showCell=True)
viewer.show()
```

## Surfaces and Slabs

### FCC Surfaces

```python
from ase.build import fcc111, fcc100, fcc110

# Au(111) surface - 4 layers, 3x3 cell, 5A vacuum
au111 = fcc111('Au', size=(3, 3, 4), vacuum=5.0)

# Pt(100) surface
pt100 = fcc100('Pt', size=(4, 4, 3), vacuum=6.0)

viewer = MolecularViewer(au111, style="metallic", showCell=True)
viewer.show()
```

### BCC Surfaces

```python
from ase.build import bcc111, bcc100, bcc110

# Fe(110) surface
fe110 = bcc110('Fe', size=(3, 3, 4), vacuum=5.0)

viewer = MolecularViewer(fe110, showCell=True)
viewer.show()
```

### General Surface

```python
from ase.build import surface

# Create any Miller index surface
atoms = bulk('Cu', 'fcc', a=3.61)
cu_211 = surface(atoms, (2, 1, 1), layers=4, vacuum=5.0)

viewer = MolecularViewer(cu_211, showCell=True)
viewer.show()
```

## Low-Dimensional Materials

### Graphene

```python
from ase.build import graphene_nanoribbon

# Zigzag nanoribbon
gnr = graphene_nanoribbon(4, 6, type='zigzag', saturated=True, vacuum=5.0)

# Armchair nanoribbon
gnr_arm = graphene_nanoribbon(4, 6, type='armchair', saturated=True, vacuum=5.0)

viewer = MolecularViewer(gnr, style="cartoon", showCell=True)
viewer.show()
```

### Carbon Nanotubes

```python
from ase.build import nanotube

# (n, m) nanotube indices
cnt_6_0 = nanotube(6, 0, length=4, vacuum=5.0)   # Zigzag
cnt_6_6 = nanotube(6, 6, length=4, vacuum=5.0)   # Armchair
cnt_8_4 = nanotube(8, 4, length=4, vacuum=5.0)   # Chiral

viewer = MolecularViewer(cnt_6_0, style="neon", backgroundColor="#000000")
viewer.show()
```

## Unit Cell Display

Toggle unit cell visibility:

```python
viewer = MolecularViewer(
    crystal,
    showCell=True,      # Show unit cell
    cellLineWidth=2.0,  # Cell line thickness
    cellColor="#888888" # Cell color
)
viewer.show()
```

## Reading Structure Files

### VASP

```python
from ase.io import read
from aseview import MolecularViewer

# Read POSCAR/CONTCAR
atoms = read("POSCAR")
viewer = MolecularViewer(atoms, showCell=True)
viewer.show()
```

```bash
aseview2 POSCAR
aseview2 CONTCAR
```

### CIF Files

```python
atoms = read("structure.cif")
viewer = MolecularViewer(atoms, showCell=True)
viewer.show()
```

```bash
aseview2 crystal.cif
```

### Other Formats

| Format | Extension | Example |
|--------|-----------|---------|
| VASP | POSCAR, CONTCAR | `aseview2 POSCAR` |
| CIF | .cif | `aseview2 structure.cif` |
| XSF | .xsf | `aseview2 charge.xsf` |
| Quantum ESPRESSO | .in | `aseview2 pw.in -f espresso-in` |
| LAMMPS | .data | `aseview2 system.data -f lammps-data` |

## Trajectory for Solid-State

### MD Trajectory

```python
from ase.io import read
from aseview import MolecularViewer

# Read VASP MD trajectory
traj = read("XDATCAR", index=":")

viewer = MolecularViewer(traj, showCell=True, showEnergyPlot=True)
viewer.show()
```

```bash
# VASP MD trajectory
aseview2 XDATCAR -i :

# ASE trajectory format
aseview2 md.traj -i :
```

### Relaxation Trajectory

```python
# Read optimization trajectory
from ase.io import read

opt_traj = read("relax.traj", index=":")

viewer = MolecularViewer(
    opt_traj,
    showCell=True,
    showEnergyPlot=True
)
viewer.show()
```

## Style Recommendations

| Structure Type | Recommended Style |
|----------------|-------------------|
| Metals | `metallic` |
| Semiconductors | `glossy` |
| Ionic crystals | `default` |
| Carbon materials | `cartoon` or `neon` |
| Surfaces | `metallic` |

## Tips for Large Systems

For systems with many atoms:

```python
viewer = MolecularViewer(
    large_system,
    atomSize=0.3,        # Smaller atoms
    bondThickness=0.08,  # Thinner bonds
    bondThreshold=0.9    # Fewer bonds detected
)
viewer.show()
```
