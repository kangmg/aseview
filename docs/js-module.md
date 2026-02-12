# ASEView.js - JavaScript Module

Standalone JavaScript library for molecular visualization. Use it in any web page without Python.

[**Live Demo**](demo.html){ .md-button }

## CDN Usage (jsDelivr)

```html
<!-- Three.js (required) -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>

<!-- ASEView.js -->
<script src="https://cdn.jsdelivr.net/gh/kangmg/aseview_v2_dev@main/aseview/static/js/styles.js"></script>
<script src="https://cdn.jsdelivr.net/gh/kangmg/aseview_v2_dev@main/aseview/static/js/aseview.js"></script>
```

## Quick Start

```html
<!DOCTYPE html>
<html>
<head>
    <style>
        #viewer { width: 600px; height: 400px; }
    </style>
</head>
<body>
    <div id="viewer"></div>

    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview_v2_dev@main/aseview/static/js/styles.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview_v2_dev@main/aseview/static/js/aseview.js"></script>

    <script>
        const viewer = new ASEView.MolecularViewer('#viewer', {
            style: 'cartoon',
            backgroundColor: '#1f2937'
        });

        // Water molecule
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

## API

### MolecularViewer

```javascript
const viewer = new ASEView.MolecularViewer(container, options);
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `style` | string | `'cartoon'` | Rendering style |
| `backgroundColor` | string | `'#1f2937'` | Background color |
| `atomSize` | number | `0.4` | Atom size scale |
| `bondThickness` | number | `0.1` | Bond thickness |
| `bondThreshold` | number | `1.2` | Bond detection threshold |
| `showBond` | boolean | `true` | Show bonds |
| `showCell` | boolean | `false` | Show unit cell |
| `rotationMode` | string | `'trackball'` | `'trackball'` or `'orbit'` |

**Styles:** `cartoon`, `glossy`, `metallic`, `neon`, `bubble`, `rowan`, `2d`, `grey`

### Methods

```javascript
// Set molecular data
viewer.setData({
    symbols: ['C', 'O', 'O'],
    positions: [[0, 0, 0], [1.16, 0, 0], [-1.16, 0, 0]],
    cell: [[10, 0, 0], [0, 10, 0], [0, 0, 10]]  // optional
});

// Update options
viewer.setOptions({ style: 'glossy' });

// Change style
viewer.setStyle('neon');

// Take screenshot
const dataUrl = viewer.screenshot();

// Clean up
viewer.dispose();
```

### NormalModeViewer

Viewer for molecular vibration animations.

```javascript
const viewer = new ASEView.NormalModeViewer('#container', {
    style: 'cartoon',
    amplitude: 1.0,
    animationSpeed: 0.05
});

viewer.setData({
    symbols: ['O', 'H', 'H'],
    positions: [[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]],
    displacements: [[0, 0, -0.07], [0, 0.42, 0.56], [0, -0.42, 0.56]]
});

// Control animation
viewer.play();
viewer.pause();
viewer.setAmplitude(2.0);
viewer.setAnimationSpeed(0.1);
```

### OverlayViewer

Compare two molecular structures by overlaying them.

```javascript
const viewer = new ASEView.OverlayViewer('#container', {
    style: 'cartoon',
    opacity1: 1.0,
    opacity2: 0.5
});

const structure1 = { symbols: ['O', 'H', 'H'], positions: [...] };
const structure2 = { symbols: ['O', 'H', 'H'], positions: [...] };

viewer.setStructures(structure1, structure2);
viewer.setOpacity(1.0, 0.3);  // Adjust opacity
```

### InteractiveViewer

Full-featured viewer with sidebar controls (same experience as Python viewer).

```javascript
const viewer = new ASEView.InteractiveViewer('#container', {
    style: 'cartoon',
    data: {
        symbols: ['O', 'H', 'H'],
        positions: [[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]]
    }
});
```

The InteractiveViewer includes:

- **Style selector**: cartoon, glossy, metallic, neon, bubble
- **Display toggles**: Show/hide bonds and unit cell
- **Size controls**: Atom size and bond thickness sliders

Also available: `InteractiveNormalModeViewer` and `InteractiveOverlayViewer` with full UI controls.

## Jekyll / GitHub Pages Usage

You can embed the viewer in Jekyll blogs (GitHub Pages):

```markdown
---
title: Molecule Example
---

# Ethanol

<div id="viewer" style="width:100%;height:400px;"></div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>
<script src="https://cdn.jsdelivr.net/gh/kangmg/aseview_v2_dev@main/aseview/static/js/styles.js"></script>
<script src="https://cdn.jsdelivr.net/gh/kangmg/aseview_v2_dev@main/aseview/static/js/aseview.js"></script>

<script>
const viewer = new ASEView.MolecularViewer('#viewer', { style: 'cartoon' });
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

Make sure `_config.yml` allows HTML:
```yaml
kramdown:
  parse_block_html: true
```
