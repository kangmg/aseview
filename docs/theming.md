# Theming

aseview supports multiple visual themes. Themes are complete HTML template sets that control the entire look and feel of the viewer — colors, typography, layout, and UI style.

## Available Themes

| Theme | Description |
|-------|-------------|
| `dark` | Default dark theme with deep grey background |
| `spring` | Light pastel theme with bright, airy colors |
| `glass` | Glassmorphism theme with frosted translucent panels |
| `darkgreen` | Dark GitHub-style theme with green accent |
| `simple` | Light macOS-style theme with colorful gradient background |

## Python API

### Per-viewer theme

Pass `theme=` to any viewer constructor:

```python
from aseview import MolecularViewer, NormalViewer, FragSelector, OverlayViewer

viewer = MolecularViewer(atoms, theme='spring')
viewer = NormalViewer(atoms, vibrations=vib, theme='spring')
viewer = FragSelector(atoms, theme='spring')
viewer = OverlayViewer([atoms1, atoms2], theme='spring')
```

### Global theme

Set once for the whole session:

```python
import aseview

aseview.set_theme('spring')         # all subsequent viewers use spring

viewer1 = MolecularViewer(atoms)    # spring
viewer2 = FragSelector(atoms)       # spring
```

### Inspect available themes

```python
aseview.list_themes()   # ['dark', 'darkgreen', 'glass', 'simple', 'spring']
aseview.get_theme()     # 'dark'  (current default)
```

## CLI

```bash
aseview molecule.xyz --theme spring
aseview molecule.xyz -t spring -o out.html
aseview molecule.xyz -v normal --theme spring --hess mol.hess
```

## JavaScript Module

### Per-viewer

```javascript
const viewer = new ASEView.MolecularViewer('#container', { theme: 'spring' });
```

### Global

```javascript
ASEView.setTheme('spring');              // all new viewers use spring
const v = new ASEView.MolecularViewer('#container');  // spring

ASEView.getTheme();                      // 'spring'
```

## VS Code Extension

If you are building a VS Code extension that wraps aseview, expose theme as a configuration setting and forward it to the viewer.

### extension.json (contributes)

```json
"contributes": {
  "configuration": {
    "title": "ASEView",
    "properties": {
      "aseview.theme": {
        "type": "string",
        "default": "dark",
        "enum": ["dark", "spring", "glass", "darkgreen", "simple"],
        "description": "Visual theme for the ASEView viewer"
      }
    }
  }
}
```

### TypeScript / JavaScript

**CLI-based extensions** — pass `--theme` when spawning aseview:

```typescript
import * as vscode from 'vscode';
import { execFile } from 'child_process';

const theme = vscode.workspace.getConfiguration('aseview').get<string>('theme', 'dark');
execFile('aseview', [filePath, '--theme', theme, '--no-browser', '-o', outPath]);
```

**Python-API-based extensions** — pass via the Python script you invoke:

```typescript
const theme = vscode.workspace.getConfiguration('aseview').get<string>('theme', 'dark');
const script = `
import aseview
viewer = aseview.MolecularViewer(atoms, theme='${theme}')
viewer.save_html('${outPath}')
`;
```

**WebView + JS module** — call `ASEView.setTheme()` before creating viewers:

```typescript
const theme = vscode.workspace.getConfiguration('aseview').get<string>('theme', 'dark');
const html = `
<script src="https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/aseview.js"></script>
<script>
  ASEView.setTheme('${theme}');
  const viewer = new ASEView.MolecularViewer('#container', { ...options });
</script>
`;
panel.webview.html = html;
```

React to configuration changes so the view updates without restarting:

```typescript
vscode.workspace.onDidChangeConfiguration(e => {
    if (e.affectsConfiguration('aseview.theme')) {
        const newTheme = vscode.workspace.getConfiguration('aseview').get<string>('theme', 'dark');
        // Re-render the webview with the new theme
        panel.webview.postMessage({ command: 'setTheme', theme: newTheme });
    }
});
```

## Adding a Custom Theme

A theme is a directory inside `aseview/themes/` containing four HTML files:

```
aseview/themes/
  my-theme/
    molecular_viewer.html
    normal_viewer.html
    overlay_viewer.html
    frag_selector.html
```

The simplest way to start is to copy an existing theme and modify the CSS variables and defaults:

```bash
cp -r aseview/themes/dark aseview/themes/my-theme
# Edit aseview/themes/my-theme/*.html
```

Once the directory exists, it will appear in `list_themes()` and be selectable via `theme='my-theme'`.
