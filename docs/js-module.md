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

        <script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
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

=== "JavaScript Only"

    ```javascript
    const viewer = new ASEView.MolecularViewer('#viewer');
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
- Energy plot visualization
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

        <script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
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

        <script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
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
- Alignment algorithms (Kabsch, Hungarian)
- RMSD calculation and display
- Multi-structure visualization

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

        <script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
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

### Data Format

```javascript
// Single structure
{
    symbols: ['C', 'H', 'O', ...],           // Required: atom symbols
    positions: [[x, y, z], ...],             // Required: atom positions (Å)
    cell: [[a1, a2, a3], [b1, b2, b3], ...], // Optional: unit cell vectors (3x3)
    energy: -123.456,                        // Optional: energy value (eV)
    forces: [[fx, fy, fz], ...],             // Optional: force vectors
    charges: [0.1, -0.2, ...],               // Optional: atomic charges
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

### Settings (Constructor Options)

All viewers accept an options object as the second argument:

```javascript
const viewer = new ASEView.MolecularViewer('#container', { style: 'neon', showBond: true });
```

Settings can also be updated at runtime:

```javascript
viewer.setSettings({ style: 'glossy', atomSize: 0.6 });
```

#### Common Settings (All Viewers)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `style` | string | `'cartoon'` | Rendering style |
| `backgroundColor` | string | `'#1f2937'` | Background color (hex) |
| `atomSize` | float | `0.4` | Atom radius scale |
| `bondThreshold` | float | `1.0` | Bond cutoff (multiple of covalent radii) |
| `bondThickness` | float | `0.1` | Bond cylinder radius |
| `showBond` | boolean | `true` | Show bonds |
| `showCell` | boolean | `true` | Show unit cell box |
| `showShading` | boolean | `true` | Enable directional lighting |
| `showAxis` | boolean | `true` | Show XYZ axis helper |
| `showForces` | boolean | `false` | Show force vectors |
| `showEnergyPlot` | boolean | `false` | Show energy vs frame plot |
| `animationSpeed` | int | `30` | Frames per second |
| `forceScale` | float | `0.5` | Force vector scale |
| `viewMode` | string | `'Perspective'` | `'Perspective'` or `'Orthographic'` |
| `rotationMode` | string | `'TrackBall'` | `'TrackBall'` or `'Orbit'` |

**Available styles:** `default`, `2d`, `cartoon`, `neon`, `glossy`, `metallic`, `rowan`, `bubble`, `grey`

#### MolecularViewer Settings

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cellColor` | string | `'#808080'` | Unit cell line color (hex) |
| `colorBy` | string | `'Element'` | `'Element'` or `'Charge'` |
| `chargeColormap` | string | `'coolwarm'` | Colormap for charge coloring |
| `normalizeCharges` | boolean | `false` | Normalize charge color range |
| `showChargeLabels` | boolean | `false` | Display charge values on atoms |

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
| `alignMolecules` | boolean | `false` | Kabsch rotation + Hungarian atom reordering |
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

        <script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
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

    <script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
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
