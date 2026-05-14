# Lightweight `view()` Helper

`view()` is the shortest Python entry point for the lightweight viewer. It is useful in notebooks when you want a PyMOL-style render surface without the full MolecularViewer control panel.

```python
from ase.build import molecule
from aseview import view

atoms = molecule("CH3CONH2")

viewer = view(
    atoms,
    styles="cinematic",
    hide_hs=False,
    center=True,
    bond_threshold=1.25,
    show_bond=True,
    high_performance=True,
    width=500,
    height=420,
)
```

<iframe src="../../assets/viewers/lite_view.html" width="100%" height="520" style="border: 1px solid #374151; border-radius: 8px;" loading="lazy"></iframe>

## Trajectory Playback

For a list of ASE `Atoms`, `view()` adds a compact play/stop button and frame slider.

```python
viewer = view(
    trajectory,
    styles="bubble",
    fps=24,
    hide_hs=True,
    center=True,
    width=500,
    height=420,
)
```

## Save the Same Viewer

`view()` returns the underlying `LiteViewer`, so the same object can be saved as standalone HTML.

```python
viewer.save_html("lite_view.html")
```

## Common Options

| Option | Description |
|--------|-------------|
| `styles` | Rendering style: `cinematic`, `bubble`, `neon`, `grey`, `2d`, `cartoon`, `glossy`, `metallic`, or `default` |
| `fps` | Playback rate for trajectories |
| `hide_hs` | Hide hydrogens without changing the underlying data |
| `center` / `centering` | Center each frame by center of mass |
| `bond_threshold` | Bond detection scale factor |
| `show_bond` | Show or hide bonds |
| `high_performance` | Use lower-overhead rendering settings for large structures |
| `width`, `height` | Notebook iframe size in pixels |
