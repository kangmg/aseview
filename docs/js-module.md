# ASEView.js - JavaScript Module

Standalone JavaScript library for molecular visualization. Use it in any web page without Python.

The JavaScript module uses the **same templates as the Python package**, ensuring identical UI and features.

[**Live Demo**](demo.html){ .md-button .md-button--primary }

## Quick Start

=== "Full HTML"

    ```html
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>ASEView Example</title>
        <style>
            #viewer { width: 100%; height: 500px; }
        </style>
    </head>
    <body>
        <div id="viewer"></div>

        <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
        <script>
            const viewer = new ASEView.MolecularViewer('#viewer', {
                viewPreset: 'top-c'
            });
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

=== "JavaScript Only"

    ```javascript
    const viewer = new ASEView.MolecularViewer('#viewer', { viewPreset: 'top-c' });
    viewer.setData({
        symbols: ['O', 'H', 'H'],
        positions: [
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469]
        ]
    });
    ```

## Available Viewers

### MolecularViewer

Full-featured molecular structure viewer with sidebar controls.

**Features:**

- Style selector (Cartoon, Glossy, Metallic, Neon, etc.)
- Bond/Cell display toggles
- Atom size and bond thickness sliders
- Trajectory animation with playback controls
- Energy and max-force plot visualization
- Radius contrast control for element-size differences
- Fixed atom highlighting from constraint metadata
- Clipboard copy as `xyz`, `extxyz`, `cif`, or `POSCAR`
- Screenshot and GIF export

=== "Full HTML"

    ```html
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>MolecularViewer Example</title>
        <style>
            #viewer { width: 100%; height: 600px; }
        </style>
    </head>
    <body>
        <div id="viewer"></div>

        <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
        <script>
            const viewer = new ASEView.MolecularViewer('#viewer');

            // Single structure (Ethanol)
            viewer.setData({
                symbols: ['C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H'],
                positions: [
                    [-0.047, 0.666, 0.000],
                    [-0.047, -0.866, 0.000],
                    [1.204, 1.144, 0.000],
                    [1.869, 0.485, 0.000],
                    [-0.570, 1.035, 0.889],
                    [-0.570, 1.035, -0.889],
                    [0.982, -1.162, 0.000],
                    [-0.557, -1.222, 0.889],
                    [-0.557, -1.222, -0.889]
                ]
            });
        </script>
    </body>
    </html>
    ```

=== "JavaScript Only"

    ```javascript
    const viewer = new ASEView.MolecularViewer('#viewer');

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

=== "Full HTML"

    ```html
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>NormalModeViewer Example</title>
        <style>
            #viewer { width: 100%; height: 600px; }
        </style>
    </head>
    <body>
        <div id="viewer"></div>

        <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
        <script>
            const viewer = new ASEView.NormalModeViewer('#viewer');

            // Equilibrium structure
            const atoms = {
                symbols: ['O', 'H', 'H'],
                positions: [
                    [0.0, 0.0, 0.117],
                    [0.0, 0.757, -0.469],
                    [0.0, -0.757, -0.469]
                ]
            };

            // Vibration data (3 modes for water)
            const vibrationData = {
                modeVectors: [
                    [[0, 0, -0.07], [0, 0.42, 0.56], [0, -0.42, 0.56]],
                    [[0, 0.07, 0], [0, -0.56, 0.42], [0, 0.56, 0.42]],
                    [[0.07, 0, 0], [-0.56, 0, 0], [-0.56, 0, 0]]
                ],
                frequencies: ["1595.32", "3657.05", "3755.93"],
                isImaginary: [false, false, false]
            };

            viewer.setVibrationData(atoms, vibrationData);
        </script>
    </body>
    </html>
    ```

=== "JavaScript Only"

    ```javascript
    const viewer = new ASEView.NormalModeViewer('#viewer');

    const atoms = {
        symbols: ['O', 'H', 'H'],
        positions: [[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]]
    };

    const vibrationData = {
        modeVectors: [
            [[0, 0, -0.07], [0, 0.42, 0.56], [0, -0.42, 0.56]],  // Mode 1
            [[0, 0.07, 0], [0, -0.56, 0.42], [0, 0.56, 0.42]],   // Mode 2
            [[0.07, 0, 0], [-0.56, 0, 0], [-0.56, 0, 0]]         // Mode 3
        ],
        frequencies: ["1595.32", "3657.05", "3755.93"],
        isImaginary: [false, false, false]
    };

    viewer.setVibrationData(atoms, vibrationData);
    ```

### OverlayViewer

Compare multiple molecular structures by overlaying them.

**Features:**

- Per-structure opacity controls
- Multi-structure visualization
- Per-structure visibility and color controls

=== "Full HTML"

    ```html
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>OverlayViewer Example</title>
        <style>
            #viewer { width: 100%; height: 600px; }
        </style>
    </head>
    <body>
        <div id="viewer"></div>

        <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
        <script>
            const viewer = new ASEView.OverlayViewer('#viewer');

            // Structure 1: Original water
            const water1 = {
                symbols: ['O', 'H', 'H'],
                positions: [
                    [0.0, 0.0, 0.117],
                    [0.0, 0.757, -0.469],
                    [0.0, -0.757, -0.469]
                ]
            };

            // Structure 2: Stretched water
            const water2 = {
                symbols: ['O', 'H', 'H'],
                positions: [
                    [0.0, 0.0, 0.050],
                    [0.0, 0.900, -0.400],
                    [0.0, -0.900, -0.400]
                ]
            };

            viewer.setStructures(water1, water2);
        </script>
    </body>
    </html>
    ```

=== "JavaScript Only"

    ```javascript
    const viewer = new ASEView.OverlayViewer('#viewer');

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

### FragSelector

Interactive fragment selector with synchronized 2D and 3D views.

**Features:**

- Click, rectangle, and lasso atom selection
- Synchronized 2D/3D fragment picking
- Copy selected or unselected atom indices
- Structure copy using the same `xyz`, `extxyz`, `cif`, and `POSCAR` formats as the main viewer

=== "Full HTML"

    ```html
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>FragSelector Example</title>
        <style>
            #viewer { width: 100%; height: 600px; }
        </style>
    </head>
    <body>
        <div id="viewer"></div>

        <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
        <script>
            const viewer = new ASEView.FragSelector('#viewer');

            viewer.setData({
                symbols: ['O', 'C', 'N', 'C', 'H', 'H', 'H', 'H', 'H'],
                positions: [
                    [ 0.424546,  1.327024,  0.008034],
                    [ 0.077158,  0.149789, -0.004249],
                    [ 0.985518, -0.878537, -0.048910],
                    [-1.371475, -0.288665, -0.000144],
                    [ 0.707952, -1.824249,  0.169942],
                    [-1.997229,  0.584922, -0.175477],
                    [-1.560842, -1.039270, -0.771686],
                    [-1.632113, -0.723007,  0.969814],
                    [ 1.953133, -0.631574,  0.111866]
                ]
            });

            viewer.setSelection([1, 2, 4, 8]);
        </script>
    </body>
    </html>
    ```

=== "JavaScript Only"

    ```javascript
    const viewer = new ASEView.FragSelector('#viewer');
    viewer.setData({
        symbols: ['O', 'C', 'N', 'C', 'H', 'H', 'H', 'H', 'H'],
        positions: [
            [ 0.424546,  1.327024,  0.008034],
            [ 0.077158,  0.149789, -0.004249],
            [ 0.985518, -0.878537, -0.048910],
            [-1.371475, -0.288665, -0.000144],
            [ 0.707952, -1.824249,  0.169942],
            [-1.997229,  0.584922, -0.175477],
            [-1.560842, -1.039270, -0.771686],
            [-1.632113, -0.723007,  0.969814],
            [ 1.953133, -0.631574,  0.111866]
        ]
    });
    viewer.setSelection([1, 2, 4, 8]);
    ```

## API Reference

### Methods

**All Viewers**

| Method | Description |
|--------|-------------|
| `setData(data)` | Set molecular data (single structure or array of structures) |
| `setSettings(options)` | Update viewer settings at runtime |
| `dispose()` | Clean up and remove viewer |

**NormalModeViewer Only**

| Method | Description |
|--------|-------------|
| `setVibrationData(atoms, vibrationData)` | Set equilibrium structure and vibration modes |

**OverlayViewer Only**

| Method | Description |
|--------|-------------|
| `setStructures(...structures)` | Set two or more structures for overlay comparison |

**FragSelector Only**

| Method | Description |
|--------|-------------|
| `setSelection(indices)` | Select atoms programmatically |
| `getSelection()` | Resolve the currently selected atom indices |
| `clearSelection()` | Clear all selected atoms |

### Data Format

```javascript
// Single structure
{
    symbols: ['C', 'H', 'O', ...],           // Required: atom symbols
    positions: [[x, y, z], ...],             // Required: atom positions (Å)
    cell: [[a1, a2, a3], [b1, b2, b3], ...], // Optional: unit cell vectors (3x3)
    pbc: [true, true, true],                  // Optional: periodic boundary flags
    energy: -123.456,                        // Optional: energy value (eV)
    forces: [[fx, fy, fz], ...],             // Optional: force vectors
    charges: [0.1, -0.2, ...],               // Optional: atomic charges
    fixed: [false, true, ...],                // Optional: fixed atom mask
    arrays: {
        move_mask: [[true, true, true], ...]  // Optional: selective dynamics mask
    },
    name: "Molecule 1"                       // Optional: display name
}

// Trajectory (array of structures)
viewer.setData([
    { symbols: [...], positions: [...], energy: -10.5 },
    { symbols: [...], positions: [...], energy: -10.3 },
    // ... more frames
]);
```

**NormalModeViewer vibration data:**

```javascript
viewer.setVibrationData(atoms, {
    modeVectors: [                 // Displacement vectors per mode [nModes][nAtoms][3]
        [[dx, dy, dz], ...],      //   Mode 0
        [[dx, dy, dz], ...],      //   Mode 1
        // ...
    ],
    frequencies: ["1595.32", "3657.05", ...],  // Frequency labels (cm⁻¹)
    isImaginary: [false, false, ...],          // Imaginary mode flags
    nFrames: 30                                // Animation frames per cycle (default: 30)
});
```

### Theming

Select a visual theme via the `theme` option or the global `setTheme()` helper:

```javascript
// Per-viewer
const viewer = new ASEView.MolecularViewer('#container', { theme: 'spring' });

// Global (affects all viewers created afterwards)
ASEView.setTheme('spring');
ASEView.getTheme();  // 'spring'

const v1 = new ASEView.MolecularViewer('#c1');    // spring
const v2 = new ASEView.NormalModeViewer('#c2');   // spring
```

| Theme | Description |
|-------|-------------|
| `dark` | Default — deep grey background |
| `spring` | Light pastel theme with bright, airy colors |
| `glass` | Glassmorphism theme with frosted translucent panels |
| `darkgreen` | Dark GitHub-style theme with green accent |
| `simple` | Light macOS-style theme with colorful gradient background |

### Settings (Constructor Options)

All viewers accept an options object as the second argument:

```javascript
const viewer = new ASEView.MolecularViewer('#container', { style: 'neon', showBond: true });
```

Settings can also be updated at runtime:

```javascript
viewer.setSettings({ style: 'glossy', atomSize: 0.6 });
```

### Camera Control

Initial camera settings can be passed in the constructor, and runtime camera
changes use `setView(viewSpec)`:

```javascript
const viewer = new ASEView.MolecularViewer('#viewer', {
    viewPreset: 'top-c',
    viewFit: 1.1
});

await viewer.setView({ preset: 'side-a', up: [0, 0, 1] });
await viewer.setView({ direction: [1, 1, 1], fit: 1.2 });
await viewer.setView({ euler: [45, 0, 30] });
await viewer.resetView();
```

`viewDirection` and preset directions are target-to-camera vectors. The camera
is placed along the normalized direction and looks back at the fitted structure
center. `viewEuler` is `[rx, ry, rz]` in degrees, XYZ order, applied to the
canonical target-to-camera vector `[0, 0, 1]`. `viewUp` is an optional up-vector
hint; invalid or parallel hints use a stable fallback. `viewFit` is a fit
multiplier.

Supported Cartesian presets are `top`, `bottom`, `front`, `back`, `left`, and
`right`. Cell-aware presets are `top-c`, `bottom-c`, `side-a`, and `side-b`.
Aliases `c` and `top` prefer the cell `c` axis when a valid cell exists, while
`a` and `b` prefer side views along the cell `a` and `b` vectors. Missing,
non-periodic, zero-length, or degenerate cells fall back to the Cartesian
direction for the requested preset.

`setView(viewSpec)` resolves with
`{ ok: true, type: 'view', view: object }`. Invalid runtime vectors reject with
an `Error` carrying `code`, `message`, `type`, and when applicable `requestId`.

### Image Export

Browser viewers can export the rendered canvas:

```javascript
const png = await viewer.savePNG({
    filename: 'structure.png',
    download: false,
    returnDataUrl: true,
    scale: 2,
    transparent: true
});

const gif = await viewer.saveGIF({
    filename: 'vibration.gif',
    frames: 30,
    delay: 40,
    sampleInterval: 20,
    download: false,
    returnDataUrl: true
});
```

Common export options are `filename`, `download` (default `true`),
`returnDataUrl` (default `false`), `scale` (default `1`), `width`, `height`,
`transparent`, and `backgroundColor`. GIF-only options are `frames`, `delay`,
and `sampleInterval`.

Successful export promises resolve with
`{ ok: true, type: 'png' | 'gif', filename, dataUrl?, width?, height? }`.
`dataUrl` is present when `returnDataUrl: true`. Failures reject with an
`Error` carrying `code`, `message`, `type`, and when applicable `requestId`.
`MolecularViewer` and `NormalModeViewer` support PNG and GIF.
`OverlayViewer` supports PNG only; `saveGIF()` rejects with
`code: 'unsupported_export'`.

Python viewers also expose headless `save_png()` and `save_gif()` when installed
with the optional export dependency:

```bash
pip install "aseview[export]"
python -m playwright install chromium
```

```python
viewer.save_png("structure.png", scale=2, transparent=False)
viewer.save_gif("trajectory.gif", frames=30, quality="high")
```

PNG quality is controlled by `scale`, `width`, and `height`. GIF
`quality="low" | "medium" | "high"` maps to the encoder sample interval.
The CLI still does not provide `--save-png` or `--save-gif`.

#### Common Settings (All Viewers)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `theme` | string | `'dark'` | Visual theme (`'dark'`, `'spring'`, `'glass'`, `'darkgreen'`, `'simple'`, …) |
| `style` | string | `'cartoon'` | Rendering style |
| `backgroundColor` | string | `'#1f2937'` | Background color (hex) |
| `atomSize` | float | `0.4` | Atom radius scale |
| `bondThreshold` | float | `1.2` | Bond cutoff (multiple of covalent radii) |
| `bondThickness` | float | `0.09` | Bond cylinder radius |
| `showBond` | boolean | `true` | Show bonds |
| `showCell` | boolean | `true` | Show unit cell box |
| `showShading` | boolean | `true` | Enable directional lighting |
| `showAxis` | boolean | `true` | Show XYZ axis helper |
| `showForces` | boolean | `false` | Show force vectors |
| `hideHydrogens` | boolean | `false` | Hide hydrogen atoms and their bonds without removing them from the data |
| `showBlur` | boolean | `false` | Apply blur to the rendered molecule canvas |
| `blurStrength` | float | `1.5` | Blur radius in pixels when `showBlur` is enabled |
| `showEnergyPlot` | boolean | `false` | Show energy vs frame plot |
| `showForceMaxPlot` | boolean | `true` | Include max force trace when forces are available |
| `radiusContrast` | float | `1.0` | Atom radius contrast (`0.0` uniform, `1.0` element radii) |
| `radiusContrastMode` | string | `'log'` | `'linear'` or `'log'` radius mapping |
| `animationSpeed` | int | `30` | Frames per second |
| `forceScale` | float | `0.5` | Force vector scale |
| `viewMode` | string | `'Perspective'` | `'Perspective'` or `'Orthographic'` |
| `viewPreset` | string \| null | `null` | Initial named camera view (`'top-c'`, `'side-a'`, `'front'`, etc.) |
| `viewDirection` | number[] \| null | `null` | Explicit target-to-camera direction vector |
| `viewEuler` | number[] \| null | `null` | XYZ Euler degrees applied to `[0, 0, 1]` |
| `viewUp` | number[] \| null | `null` | Optional camera up-vector hint |
| `viewFit` | number | `1` | Camera fit multiplier |
| `rotationMode` | string | `'TrackBall'` | `'TrackBall'` or `'Orbit'` |

**Available styles:** `default`, `2d`, `cartoon`, `neon`, `glossy`, `metallic`, `cinematic`, `rowan`, `bubble`, `grey`

#### MolecularViewer Settings

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cellColor` | string | `'#808080'` | Unit cell line color (hex) |
| `colorBy` | string | `'Element'` | `'Element'` or `'Charge'` |
| `chargeColormap` | string | `'coolwarm'` | Colormap for charge coloring |
| `normalizeCharges` | boolean | `false` | Normalize charge color range |
| `showChargeLabels` | boolean | `false` | Display charge values on atoms |
| `showConstraint` | boolean | `false` | Highlight fixed atoms with a yellow overlay |

#### NormalModeViewer Settings

Inherits all common settings, plus:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `showModeVector` | boolean | `false` | Display vibration displacement arrows |
| `showShadow` | boolean | `false` | Show shadow beneath molecule |

Vibration playback parameters (set via `setVibrationData`):

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `amplitude` | float | `0.75` | Displacement amplitude scale |
| `nFrames` | int | `30` | Animation frames per cycle |

#### OverlayViewer Settings

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `colorBy` | string | `'Atom'` | `'Atom'`, `'Molecule'`, or `'Colormap'` |
| `colormap` | string | `'viridis'` | Colormap name (when `colorBy: 'Colormap'`) |
| `showShading` | boolean | `false` | Shading (default differs from other viewers) |

**Available colormaps:** `viridis`, `plasma`, `coolwarm`, `jet`, `rainbow`, `grayscale`

## Jekyll / GitHub Pages

=== "Full HTML"

    ```html
    ---
    title: My Molecule
    ---
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            #viewer { width: 100%; height: 500px; }
        </style>
    </head>
    <body>
        <div id="viewer"></div>

        <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
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
    </body>
    </html>
    ```

=== "Embed in Markdown"

    ```html
    ---
    title: My Molecule
    ---

    # My Research

    Some text about the molecule...

    <div id="viewer" style="width:100%; height:500px;"></div>

    <script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
    <script>
    const viewer = new ASEView.MolecularViewer('#viewer');
    viewer.setData({
        symbols: ['O', 'H', 'H'],
        positions: [[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]]
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
