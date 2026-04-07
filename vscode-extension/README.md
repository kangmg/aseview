# ASEView VS Code Extension (Read-only Viewer)

This extension opens ASE-supported structure files in a VS Code native custom editor tab.

- Parsing: `ase-ts`
- Rendering UI: `ASEView.js` templates
- No Python runtime required

## What it does

- Opens structures in a VS Code webview tab (not an external browser).
- Supports Explorer context menu: `Open With ASEView`.
- Supports command palette:
  - `ASEView: Open Active File`
  - `ASEView: Refresh Current View`

## Supported default file patterns

- `*.xyz`
- `*.extxyz`
- `*.cif`
- `*.pdb`
- `*.vasp`
- `POSCAR`
- `CONTCAR`

## Build and package

```bash
cd vscode-extension
npm install
npm run compile
npm run package
```

Install the resulting `.vsix` in VS Code:

- Extensions view (`Ctrl+Shift+X`)
- `...` menu
- `Install from VSIX...`

## Settings

- `aseview.defaultViewer`: `auto | molecular | overlay | frag | normal`
- `aseview.defaultStyle`: viewer style
- `aseview.readIndex`: ASE-style index/slice (default `:`)
- `aseview.readFormat`: optional parser format override
