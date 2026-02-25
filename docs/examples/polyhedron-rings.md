# Polyhedron & Ring Highlight

Two visualization features for structural analysis:

- **Polyhedron** — render coordination polyhedra around atoms in solid-state structures
- **Ring Highlight** — fill aromatic and aliphatic rings in molecules with a semi-transparent color

---

## Polyhedron Visualization

### Live Demo: CaSiO₃ (Tetragonal Perovskite)

SiO₄ tetrahedra and CaO₈ polyhedra in CaSiO₃ (2×2×2 supercell, space group I4/mcm):

<iframe src="../../assets/viewers/casio3_polyhedron.html" width="100%" height="520" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

### How It Works

For each atom with coordination number (CN) **≥ 4**, a convex hull polyhedron is drawn around its bonded neighbors. The face color is the **average CPK color** of the three corner atoms for each triangular face.

The bond detection uses **scaled covalent radii** (controlled by `bondThreshold`); the polyhedra follow exactly the same neighbor list.

### Code Example

=== "From CIF"

    ```python
    from ase.io import read
    from aseview import MolecularViewer

    # Read CaSiO3 crystal structure from CIF
    atoms = read("CaSiO3.cif")
    atoms = atoms * (2, 2, 2)  # supercell

    viewer = MolecularViewer(
        atoms,
        style="default",
        showCell=True,
        showPolyhedron=True,       # enable polyhedra
        polyhedronOpacity=0.25,    # 0.0–1.0
        bondThreshold=1.1,
    )
    viewer.show()
    ```

=== "Inline"

    ```python
    from ase.build import bulk
    from aseview import MolecularViewer

    # NaCl rocksalt — Na is octahedrally coordinated (CN=6)
    atoms = bulk("NaCl", "rocksalt", a=5.64) * (2, 2, 2)

    viewer = MolecularViewer(
        atoms,
        showPolyhedron=True,
        polyhedronOpacity=0.3,
    )
    viewer.show()
    ```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `showPolyhedron` | `bool` | `False` | Enable polyhedron rendering |
| `polyhedronOpacity` | `float` | `0.25` | Face transparency (0 = invisible, 1 = solid) |

!!! note "Coordination threshold"
    Only atoms with **CN ≥ 4** generate polyhedra. Atoms with fewer bonds (e.g., terminal O in a chain) are skipped.

### UI Controls

In the sidebar → **Display Settings**:

- **Polyhedron** toggle — show/hide all polyhedra
- **Polyhedron Opacity** slider — adjust transparency in real time

---

## Ring Highlight

### Live Demo: Carbazole

Tricyclic aromatic system: two benzene rings fused to a pyrrole ring (C₁₂H₉N core):

<iframe src="../../assets/viewers/carbazole_rings.html" width="100%" height="520" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

### How It Works

A BFS-based ring detection algorithm finds all **smallest rings of 4–8 atoms** in the molecular bond graph. For each ring, a fan-triangulated face is drawn using the CPK color of the **most common element** in that ring.

### Code Example

=== "From XYZ"

    ```python
    from ase.io import read
    from aseview import MolecularViewer

    # Carbazole from XYZ file
    atoms = read("carbazole.xyz")

    viewer = MolecularViewer(
        atoms,
        style="cartoon",
        showRings=True,       # enable ring highlighting
        ringOpacity=0.35,     # 0.0–1.0
        bondThreshold=1.15,
    )
    viewer.show()
    ```

=== "Inline (benzene)"

    ```python
    from ase.build import molecule
    from aseview import MolecularViewer

    benzene = molecule("C6H6")

    viewer = MolecularViewer(
        benzene,
        showRings=True,
        ringOpacity=0.4,
    )
    viewer.show()
    ```

=== "Naphthalene"

    ```python
    from ase.build import molecule
    from aseview import MolecularViewer

    naphthalene = molecule("C10H8")

    viewer = MolecularViewer(
        naphthalene,
        style="cartoon",
        showRings=True,
        ringOpacity=0.4,
    )
    viewer.show()
    ```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `showRings` | `bool` | `False` | Enable ring face rendering |
| `ringOpacity` | `float` | `0.30` | Face transparency (0 = invisible, 1 = solid) |

!!! note "Ring size range"
    Rings of **4–8 atoms** are detected. Larger macrocycles (crown ethers, porphyrins > 8 atoms) are not highlighted.

### UI Controls

In the sidebar → **Display Settings**:

- **Ring Highlight** toggle — show/hide ring faces
- **Ring Opacity** slider — adjust transparency in real time

---

## Combining Both Features

Both features can be enabled simultaneously:

```python
from ase.io import read
from aseview import MolecularViewer

# Any structure with rings and coordination environments
atoms = read("your_structure.cif")

viewer = MolecularViewer(
    atoms,
    showPolyhedron=True,
    polyhedronOpacity=0.2,
    showRings=True,
    ringOpacity=0.35,
    showCell=True,
)
viewer.show()
```
