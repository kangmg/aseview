# Basic Visualization

This example demonstrates basic molecular visualization with different styles.

## Water Molecule

=== "Python"

    ```python
    from ase.build import molecule
    from aseview import MolecularViewer

    water = molecule('H2O')
    viewer = MolecularViewer(water)
    viewer.show()
    ```

=== "CLI"

    ```bash
    aseview2 water.xyz
    ```

<div class="viewer-container" id="water-viewer">
    <iframe
        src="../assets/viewers/water.html"
        width="100%"
        height="500px"
        style="border: 1px solid #374151; border-radius: 8px;"
        loading="lazy">
    </iframe>
</div>

## Ethanol with Different Styles

### Cartoon Style (Default)

```python
viewer = MolecularViewer(atoms, style="cartoon")
```

### Neon Style

```python
viewer = MolecularViewer(atoms, style="neon")
```

### Glossy Style

```python
viewer = MolecularViewer(atoms, style="glossy")
```

## Customizing Appearance

```python
viewer = MolecularViewer(
    atoms,
    atomSize=0.5,           # Larger atoms
    bondThickness=0.15,     # Thicker bonds
    backgroundColor="#000000",  # Black background
    style="metallic"
)
viewer.show()
```

## Saving to HTML

```python
viewer = MolecularViewer(atoms)
viewer.save_html("molecule.html")
```

This creates a standalone HTML file that can be opened in any browser without requiring Python.
