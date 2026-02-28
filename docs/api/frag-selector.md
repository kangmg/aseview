# FragSelector

Interactive fragment selector with synchronized 2D (SVG) and 3D (THREE.js) views.
Click any atom in either panel to toggle its selection.
Selected and unselected atom indices are displayed live in the sidebar and can be copied to the clipboard with one click.

## Constructor

```python
FragSelector(data, **kwargs)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `data` | `Atoms` or `str` | A single ASE Atoms object or a path to a structure file. Multi-frame files are accepted — only the **first frame** is used. | Required |

### Keyword Arguments (Settings)

#### Geometry

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `bondThreshold` | `float` | Scale factor applied to the sum of covalent radii for bond detection | `1.2` |
| `atomSize` | `float` | Relative atom sphere size in the 3D view | `0.4` |

#### Style

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `style` | `str` | 3D rendering style (see table below) | `"cartoon"` |
| `backgroundColor` | `str` | Background colour of the 3D panel (hex string) | `"#1f2937"` |

##### Available 3D Styles

| Value | Description |
|-------|-------------|
| `"cartoon"` | Cartoon style with black bonds (default) |
| `"default"` | Standard CPK coloring |
| `"glossy"` | Shiny reflective surface |
| `"metallic"` | Metallic appearance |
| `"bubble"` | Bubble-like appearance |
| `"rowan"` | Rowan-inspired style |
| `"neon"` | Glowing neon effect |

#### Rendering

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `showBond` | `bool` | Draw bonds in the 3D view | `True` |
| `showShading` | `bool` | Enable directional shading on atom spheres | `True` |

## Methods

### show()

Display the viewer inline in a Jupyter notebook.

```python
viewer.show(width='100%', height=600)
```

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `width` | `str` or `int` | Width of the iframe (`"100%"`, `800`, …) | `"100%"` |
| `height` | `int` | Height of the iframe in pixels | `600` |

### get_html()

Return the fully self-contained HTML as a string.
All JavaScript dependencies are inlined — no internet connection required to open the file.

```python
html_string = viewer.get_html()
```

### save_html()

Save the viewer as a standalone HTML file.

```python
viewer.save_html("selector.html")
```

## UI Overview

The viewer is divided into three areas:

```
┌──────────┬──────────────────┬───┬──────────────────┐
│ Sidebar  │   2D Structure   │ ║ │   3D Structure   │
│          │   (SVG panel)    │   │  (THREE.js panel) │
│Selection │                  │   │                   │
│Display   │  click to select │   │  click to select  │
└──────────┴──────────────────┴───┴──────────────────┘
```

### 2D Panel (SVG)

| Action | Effect |
|--------|--------|
| **Click** atom | Toggle selection |
| **Drag** background | Pan view |
| **Scroll** | Zoom in / out |
| **Double-click** background | Reset fit |
| **Hover** atom | Show element + index tooltip |

The 2D layout is computed automatically from the 3D coordinates using a **PCA → Fruchterman-Reingold** pipeline:

1. Bonds detected with the same covalent-radii threshold as the 3D view
2. Connected components separated by BFS
3. Each component projected to 2D via PCA, then refined with Fruchterman-Reingold force-directed layout
4. Average bond length normalised to 1.5 layout units
5. Fragments arranged left-to-right, largest fragment first

Bond lines use split colours (each half inherits the endpoint atom's CPK colour), matching the appearance of the 3D view.

### 3D Panel (THREE.js)

| Action | Effect |
|--------|--------|
| **Click** atom | Toggle selection |
| **Left-drag** | Rotate (TrackballControls) |
| **Right-drag** | Pan |
| **Scroll** | Zoom |

Selected atoms are highlighted with a semi-transparent blue sphere overlay (no re-render of the full molecule required).

### Splitter

Drag the vertical divider between the two panels to resize them.

### Sidebar — Selection card

| Element | Description |
|---------|-------------|
| **Selected** counter + textarea | Live count and `[i, j, ...]` list of selected indices |
| **Unselected** counter + textarea | Live count and `[i, j, ...]` list of unselected indices |
| **Copy** | Copy the corresponding index list to the clipboard |
| **Select All** | Add all atoms to the selection |
| **Clear** | Remove all atoms from the selection |
| **Invert** | Swap selected ↔ unselected |

### Sidebar — Display card

| Control | Description |
|---------|-------------|
| Bond Threshold slider | Adjust bond detection cutoff (triggers layout + 3D re-render) |
| Atom Size slider | Scale atom sphere radius in the 3D view |
| 3D Style selector | Switch between visual styles |
| Background colour picker | Change 3D panel background |
