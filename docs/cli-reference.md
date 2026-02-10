# CLI Reference

## Synopsis

```bash
aseview2 [OPTIONS] FILES...
```

## Arguments

| Argument | Description |
|----------|-------------|
| `FILES` | One or more input files in any ASE-supported format (xyz, cif, pdb, vasp, etc.) |

## Options

### Input Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--index` | `-i` | Index or slice for trajectory selection | All frames |
| `--format` | `-f` | Force file format (auto-detected if not specified) | Auto |

#### Index Examples

| Value | Description |
|-------|-------------|
| `0` | First frame only |
| `-1` | Last frame only |
| `:` | All frames |
| `0:10` | Frames 0-9 |
| `::2` | Every 2nd frame |
| `10:20:2` | Frames 10-19, step 2 |
| `-10:` | Last 10 frames |

### Viewer Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--viewer` | `-v` | Viewer type: `molecular`, `overlay`, `normal` | Auto-selected |

#### Viewer Auto-Selection

- Multiple files → `overlay`
- Single file with `--hess` → `normal`
- Otherwise → `molecular`

### Style Options

| Option | Description | Default |
|--------|-------------|---------|
| `--style` | Visual style | `cartoon` |

#### Available Styles

| Style | Description |
|-------|-------------|
| `default` | Standard CPK coloring |
| `cartoon` | Cartoon style with black bonds |
| `neon` | Glowing neon effect |
| `glossy` | Shiny reflective surface |
| `metallic` | Metallic appearance |
| `rowan` | Rowan-inspired style |
| `grey` | Greyscale rendering |

### Overlay Options

| Option | Description | Default |
|--------|-------------|---------|
| `--cmap` | Colormap for gradient coloring | None |

#### Available Colormaps

| Colormap | Description |
|----------|-------------|
| `viridis` | Perceptually uniform (blue → green → yellow) |
| `plasma` | Perceptually uniform (purple → orange → yellow) |
| `coolwarm` | Diverging (blue → white → red) |
| `jet` | Rainbow (blue → cyan → yellow → red) |
| `rainbow` | Full spectrum |
| `grayscale` | Black to white |

### Normal Mode Options

| Option | Description | Default |
|--------|-------------|---------|
| `--hess` | Path to Hessian file (ORCA .hess format) | None |

### Output Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--output` | `-o` | Save HTML to file instead of serving | None |
| `--port` | `-p` | HTTP server port | 8080 |
| `--no-browser` | | Don't open browser automatically | False |

### Utility Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--kill` | `-k` | Kill existing process on port before starting | False |
| `--help` | `-h` | Show help message | |

## Examples

### Basic Viewing

```bash
# Single structure
aseview2 molecule.xyz

# Crystal structure
aseview2 structure.cif

# Specific style
aseview2 molecule.xyz --style neon
```

### Trajectory

```bash
# All frames with animation
aseview2 trajectory.xyz

# Subset of frames
aseview2 trajectory.xyz -i 0:100:5

# Save as HTML
aseview2 trajectory.xyz -o trajectory.html
```

### Overlay Comparison

```bash
# Compare two structures
aseview2 reactant.xyz product.xyz

# Trajectory overlay with colormap
aseview2 optimization.xyz -v overlay --cmap viridis

# Multiple files
aseview2 conf1.xyz conf2.xyz conf3.xyz --cmap plasma
```

### Normal Modes

```bash
# ORCA Hessian
aseview2 molecule.xyz --hess orca.hess

# Explicit viewer type
aseview2 molecule.xyz --hess orca.hess -v normal
```

### Server Options

```bash
# Custom port
aseview2 molecule.xyz -p 9000

# Kill existing server first
aseview2 molecule.xyz -k -p 8080

# No auto-open browser
aseview2 molecule.xyz --no-browser
```

### SSH Forwarding

```bash
# Server side
aseview2 molecule.xyz -p 8080 --no-browser

# Client side
ssh -L 8080:localhost:8080 user@server
# Then open http://localhost:8080 in browser
```
