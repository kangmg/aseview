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
