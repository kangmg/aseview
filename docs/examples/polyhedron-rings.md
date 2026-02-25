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

Optionally, thin **edge lines** (`THREE.EdgesGeometry`) can be overlaid on each face to make the geometry more legible.

### Code Example

=== "NaCl (bulk)"

    ```python
    from ase.build import bulk
    from aseview import MolecularViewer

    # NaCl rocksalt — Na & Cl are octahedrally coordinated (CN=6)
    atoms = bulk("NaCl", "rocksalt", a=5.64) * (2, 2, 2)

    viewer = MolecularViewer(
        atoms,
        showPolyhedron=True,       # enable polyhedra
        polyhedronOpacity=0.25,    # face transparency 0.0–1.0
        showPolyhedronEdge=True,   # dark wireframe edges
        polyhedronEdgeOpacity=0.7, # edge visibility 0.0–1.0
    )
    viewer.show()
    ```

=== "From CIF"

    ```python
    from ase.io import read
    from aseview import MolecularViewer

    atoms = read("structure.cif") * (2, 2, 2)

    viewer = MolecularViewer(
        atoms,
        style="default",
        showCell=True,
        showPolyhedron=True,
        polyhedronOpacity=0.25,
        bondThreshold=1.1,
    )
    viewer.show()
    ```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `showPolyhedron` | `bool` | `False` | Enable polyhedron rendering |
| `polyhedronOpacity` | `float` | `0.25` | Face transparency (0 = invisible, 1 = solid) |
| `showPolyhedronEdge` | `bool` | `True` | Show wireframe edge lines on each face |
| `polyhedronEdgeOpacity` | `float` | `0.70` | Edge line opacity (0 = invisible, 1 = opaque) |

!!! note "Coordination threshold"
    Only atoms with **CN ≥ 4** generate polyhedra. Atoms with fewer bonds (e.g., terminal O in a chain) are skipped.

!!! note "Edge line width"
    WebGL does not support `linewidth > 1`. Edges are always 1 px wide; adjust `polyhedronEdgeOpacity` to make them more or less prominent.

### UI Controls

When **Polyhedron** is toggled **ON** in **Display Settings**, the following controls appear:

| Control | Description |
|---------|-------------|
| **Polyhedron Opacity** slider | Adjust face transparency |
| **Polyhedron Edge** toggle | Show / hide wireframe edges |
| **Edge Opacity** slider | Adjust edge line visibility |
| **Exclude Elements** chips | Click a chip to remove that element's polyhedra from the view; click again to restore |

---

## Ring Highlight

### Live Demo: Carbazole

Tricyclic aromatic system: two benzene rings fused to a pyrrole ring (C₁₂H₉N core):

<iframe src="../../assets/viewers/carbazole_rings.html" width="100%" height="520" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

---

### How It Works

A BFS-based ring detection algorithm finds all **smallest rings of 4–8 atoms** in the molecular bond graph. For each ring, a fan-triangulated face is drawn using the CPK color of the **most common element** in that ring.

### Code Example

=== "Benzene"

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

=== "From XYZ"

    ```python
    from ase.io import read
    from aseview import MolecularViewer

    atoms = read("carbazole.xyz")

    viewer = MolecularViewer(
        atoms,
        style="cartoon",
        showRings=True,
        ringOpacity=0.35,
        bondThreshold=1.15,
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

atoms = read("your_structure.cif")

viewer = MolecularViewer(
    atoms,
    showPolyhedron=True,
    polyhedronOpacity=0.2,
    showPolyhedronEdge=True,
    polyhedronEdgeOpacity=0.7,
    showRings=True,
    ringOpacity=0.35,
    showCell=True,
)
viewer.show()
```
