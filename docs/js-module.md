# ASEView.js - JavaScript Module

Standalone JavaScript library for molecular visualization. Use it in any web page without Python.

[**Live Demo**](demo.html){ .md-button }

## CDN Usage (GitHub Raw)

```html
<!-- Three.js (required) -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>

<!-- ASEView.js -->
<script src="https://raw.githubusercontent.com/kangmg/aseview_v2_dev/main/test_js_module/dist/aseview.umd.min.js"></script>
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
    <script src="https://raw.githubusercontent.com/kangmg/aseview_v2_dev/main/test_js_module/dist/aseview.umd.min.js"></script>

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
| `showBond` | boolean | `true` | Show bonds |
| `showCell` | boolean | `false` | Show unit cell |
| `atomSize` | number | `1.0` | Atom size scale |
| `bondThickness` | number | `0.15` | Bond thickness |
| `bondThreshold` | number | `1.2` | Bond detection threshold |
| `colorBy` | string | `'Element'` | Color mode: `'Element'` or `'Charge'` |

**Styles:** `cartoon`, `glossy`, `metallic`, `neon`, `bubble`, `rowan`, `2d`, `grey`

### Methods

```javascript
// Set molecular data
viewer.setData({
    symbols: ['C', 'O', 'O'],
    positions: [[0, 0, 0], [1.16, 0, 0], [-1.16, 0, 0]],
    charges: [0.4, -0.2, -0.2],  // optional
    cell: [[10, 0, 0], [0, 10, 0], [0, 0, 10]]  // optional
});

// Update options
viewer.setOptions({ style: 'glossy', colorBy: 'Charge' });

// Take screenshot
const dataUrl = viewer.screenshot();

// Clean up
viewer.dispose();
```

## Example: Caffeine with Charge Coloring

```html
<!DOCTYPE html>
<html>
<head>
    <style>
        body { margin: 0; background: #111; }
        #viewer { width: 100vw; height: 100vh; }
        .controls {
            position: fixed;
            top: 20px;
            left: 20px;
            z-index: 100;
        }
        button {
            background: #333;
            color: #fff;
            border: none;
            padding: 8px 16px;
            margin: 4px;
            cursor: pointer;
            border-radius: 4px;
        }
        button:hover { background: #555; }
        button.active { background: #3b82f6; }
    </style>
</head>
<body>
    <div class="controls">
        <button onclick="setStyle('cartoon')" class="active">Cartoon</button>
        <button onclick="setStyle('glossy')">Glossy</button>
        <button onclick="setStyle('neon')">Neon</button>
        <button onclick="viewer.setOptions({colorBy: 'Charge'})">Charge</button>
    </div>
    <div id="viewer"></div>

    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>
    <script src="https://raw.githubusercontent.com/kangmg/aseview_v2_dev/main/test_js_module/dist/aseview.umd.min.js"></script>

    <script>
        const viewer = new ASEView.MolecularViewer('#viewer', {
            style: 'cartoon',
            backgroundColor: '#111827'
        });

        // Caffeine molecule
        viewer.setData({
            symbols: ['N', 'C', 'N', 'C', 'C', 'N', 'C', 'N', 'C', 'C', 'O', 'N', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
            positions: [
                [1.287, -1.067, 0.000], [1.998, 0.121, 0.000], [1.282, 1.294, 0.000],
                [-0.090, 1.050, 0.000], [-0.462, -0.326, 0.000], [0.133, -1.320, 0.000],
                [-0.678, 2.247, 0.000], [3.414, 0.172, 0.000], [-1.771, -0.809, 0.000],
                [-2.781, 0.187, 0.000], [-0.114, 3.308, 0.000], [-2.009, 2.007, 0.000],
                [-2.549, 1.536, 0.000], [-2.060, -2.000, 0.000], [3.720, 1.145, 0.000],
                [3.774, -0.332, 0.890], [3.774, -0.332, -0.890], [-2.732, 2.828, 0.000],
                [-0.393, -2.247, 0.000], [-3.598, 1.832, 0.000], [1.769, -1.964, 0.000],
                [-2.883, -0.101, 1.029], [-2.883, -0.101, -1.029], [-3.803, -0.166, 0.000]
            ],
            charges: [
                -0.360, 0.520, -0.560, 0.320, 0.180, -0.280,
                0.540, -0.120, 0.220, 0.180, -0.420, -0.340,
                0.280, -0.380, 0.080, 0.060, 0.060, 0.160,
                0.120, 0.060, 0.220, 0.040, 0.040, 0.040
            ]
        });

        function setStyle(style) {
            viewer.setOptions({ style });
            document.querySelectorAll('button').forEach(btn => {
                if (['cartoon', 'glossy', 'neon'].includes(btn.textContent.toLowerCase())) {
                    btn.classList.toggle('active', btn.textContent.toLowerCase() === style);
                }
            });
        }
    </script>
</body>
</html>
```

## Local Development

```bash
cd test_js_module
npm install
npm run build
```

Build outputs:
- `dist/aseview.umd.js` - UMD bundle (browsers)
- `dist/aseview.umd.min.js` - Minified UMD
- `dist/aseview.esm.js` - ES module
