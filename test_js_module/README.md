# ASEView.js

JavaScript library for molecular visualization, designed to be embedded via CDN.

## Quick Start

```html
<!-- Include Three.js -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>

<!-- Include ASEView.js -->
<script src="https://cdn.jsdelivr.net/npm/aseview-js/dist/aseview.umd.min.js"></script>

<div id="viewer" style="width: 600px; height: 400px;"></div>

<script>
  const viewer = new ASEView.MolecularViewer('#viewer', {
    style: 'cartoon',
    colorBy: 'Element'
  });

  viewer.setData({
    symbols: ['O', 'H', 'H'],
    positions: [
      [0, 0, 0],
      [0.96, 0, 0],
      [-0.24, 0.93, 0]
    ]
  });
</script>
```

## Development

```bash
# Install dependencies
npm install

# Build
npm run build

# Serve examples
npm run serve
# Then open http://localhost:8080/examples/demo.html
```

## API

### MolecularViewer

```javascript
const viewer = new ASEView.MolecularViewer(container, options);
```

**Options:**
- `style` - Rendering style: `'cartoon'`, `'glossy'`, `'metallic'`, `'neon'`, `'bubble'`, `'2d'`
- `backgroundColor` - Background color (default: `'#1f2937'`)
- `colorBy` - Color mode: `'Element'` or `'Charge'`
- `atomSize` - Atom size scale (default: `1.0`)
- `bondThickness` - Bond radius (default: `0.15`)
- `bondThreshold` - Bond detection threshold (default: `1.5`)
- `showBond` - Show bonds (default: `true`)
- `showCell` - Show unit cell (default: `false`)

**Methods:**
- `setData(data)` - Set molecular data
- `setOptions(options)` - Update viewer options
- `screenshot()` - Get PNG data URL
- `dispose()` - Clean up resources

**Data format:**
```javascript
{
  symbols: ['C', 'H', 'H', 'H', 'H'],  // Element symbols
  positions: [[0,0,0], [1,0,0], ...],   // [x, y, z] coordinates
  charges: [0.1, -0.02, ...],           // Optional: partial charges
  cell: [[a1,a2,a3], [b1,b2,b3], ...]   // Optional: unit cell vectors
}
```

## File Structure

```
test_js_module/
├── src/
│   ├── core/
│   │   └── Scene.js          # Three.js scene management
│   ├── viewers/
│   │   └── MolecularViewer.js
│   ├── renderers/
│   │   ├── AtomRenderer.js
│   │   └── BondRenderer.js
│   ├── utils/
│   │   └── atomInfo.js       # Atom colors, radii, bond detection
│   └── index.js              # Entry point
├── dist/                      # Built files
├── examples/
│   └── demo.html
├── package.json
└── rollup.config.js
```
