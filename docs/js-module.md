# ASEView.js - JavaScript Module

Standalone JavaScript library for molecular visualization. Use it in any web page without Python.

The JavaScript module uses the **same templates as the Python package**, ensuring identical UI and features.

[**Live Demo**](demo.html){ .md-button }

## Quick Start

```html
<!DOCTYPE html>
<html>
<head>
    <style>
        #viewer { width: 100%; height: 500px; }
    </style>
</head>
<body>
    <div id="viewer"></div>

    <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview_v2_dev@main/aseview/static/js/aseview.js"></script>
    <script>
        const viewer = new ASEView.MolecularViewer('#viewer');
        viewer.setData({
            symbols: ['O', 'H', 'H'],
            positions: [
                [0.0, 0.0, 0.117],
                [0.0, 0.757, -0.469],
                [0.0, -0.757, -0.469]
            ]
        });
    </script>
</body>
</html>
```

## Available Viewers

### MolecularViewer

Full-featured molecular structure viewer with sidebar controls.

**Features:**

- Style selector (Cartoon, Glossy, Metallic, Neon, etc.)
- Bond/Cell display toggles
- Atom size and bond thickness sliders
- Trajectory animation with playback controls
- Energy plot visualization
- Screenshot and GIF export

```javascript
const viewer = new ASEView.MolecularViewer('#container');

// Single structure
viewer.setData({
    symbols: ['C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H'],
    positions: [
        [-0.047, 0.666, 0.000],
        [-0.047, -0.866, 0.000],
        [1.204, 1.144, 0.000],
        // ... more positions
    ]
});

// Trajectory (array of structures)
viewer.setData([
    { symbols: [...], positions: [...], energy: -10.5 },
    { symbols: [...], positions: [...], energy: -10.3 },
    // ... more frames
]);
```

### NormalModeViewer

Viewer for molecular vibration animations.

**Features:**

- Mode selector dropdown
- Amplitude control slider
- Animation speed control
- Play/Pause toggle
- Frequency display

```javascript
const viewer = new ASEView.NormalModeViewer('#container');

// With vibration data
viewer.setVibrationData(
    // Equilibrium structure
    {
        symbols: ['O', 'H', 'H'],
        positions: [[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]]
    },
    // Vibration data
    {
        modeVectors: [
            [[0, 0, -0.07], [0, 0.42, 0.56], [0, -0.42, 0.56]],  // Mode 1
            [[0, 0.58, -0.39], [0, -0.42, 0.56], [0, -0.42, 0.56]],  // Mode 2
            // ... more modes
        ],
        frequencies: ["1595.32", "3657.05", "3755.93"],
        isImaginary: [false, false, false]
    }
);
```

### OverlayViewer

Compare multiple molecular structures by overlaying them.

**Features:**

- Per-structure opacity controls
- Alignment algorithms (Kabsch, Hungarian)
- RMSD calculation and display
- Multi-structure visualization

```javascript
const viewer = new ASEView.OverlayViewer('#container');

const structure1 = {
    symbols: ['O', 'H', 'H'],
    positions: [[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]]
};

const structure2 = {
    symbols: ['O', 'H', 'H'],
    positions: [[0, 0, 0.050], [0, 0.900, -0.400], [0, -0.900, -0.400]]
};

viewer.setStructures(structure1, structure2);

// Or use setData with an array
viewer.setData([structure1, structure2]);
```

## API Reference

### Constructor Options

All viewers accept an options object:

```javascript
const viewer = new ASEView.MolecularViewer('#container', {
    style: 'Cartoon',           // Rendering style
    backgroundColor: '#1f2937', // Background color
    atomSize: 0.4,              // Atom size scale
    bondThickness: 0.1,         // Bond thickness
    showBond: true,             // Show bonds
    showCell: false             // Show unit cell
});
```

### Methods

| Method | Description |
|--------|-------------|
| `setData(data)` | Set molecular data (structure or trajectory) |
| `setSettings(options)` | Update viewer settings |
| `dispose()` | Clean up and remove viewer |

### Data Format

```javascript
{
    symbols: ['C', 'H', 'O', ...],           // Required: atom symbols
    positions: [[x, y, z], ...],             // Required: atom positions (Å)
    cell: [[a1, a2, a3], [b1, b2, b3], ...], // Optional: unit cell vectors
    energy: -123.456,                        // Optional: energy value (eV)
    forces: [[fx, fy, fz], ...],             // Optional: force vectors
    charges: [0.1, -0.2, ...]                // Optional: atomic charges
}
```

## Jekyll / GitHub Pages

Embed viewers in Jekyll blogs or GitHub Pages:

```html
---
title: My Molecule
---

<div id="viewer" style="width:100%; height:500px;"></div>

<script src="https://cdn.jsdelivr.net/gh/kangmg/aseview_v2_dev@main/aseview/static/js/aseview.js"></script>
<script>
const viewer = new ASEView.MolecularViewer('#viewer');
viewer.setData({
    symbols: ['C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H'],
    positions: [
        [-0.047, 0.666, 0.000], [-0.047, -0.866, 0.000],
        [1.204, 1.144, 0.000], [1.869, 0.485, 0.000],
        [-0.570, 1.035, 0.889], [-0.570, 1.035, -0.889],
        [0.982, -1.162, 0.000], [-0.557, -1.222, 0.889],
        [-0.557, -1.222, -0.889]
    ]
});
</script>
```

## Architecture

The JavaScript module uses the **same HTML templates** as the Python package:

```
aseview.js (thin wrapper)
    ↓
iframe loads template from CDN
    ↓
postMessage sends data to template
    ↓
Template renders (same code as Python viewer)
```

This ensures:

- **Consistent UI**: Web and Python viewers look identical
- **Single source of truth**: One template, no duplicate code
- **Automatic updates**: Python template improvements apply to web
