# VS CODE EXTENSION GUIDE

## OVERVIEW
This directory is a VS Code custom editor extension for molecular structure
files. It parses and writes structures with Node workers around `ase-ts`, then
embeds the aseview browser viewer in a VS Code webview.
The user-facing extension is currently documented as a read-only viewer.

## STRUCTURE
```
vscode-extension/
  src/extension.ts          # Custom editor provider and webview bridge
  node/parse_with_ase_ts.mjs # Structure file to viewer-frame JSON
  node/write_with_ase_ts.mjs # Viewer-frame JSON to copy/export formats
  media/                    # Bundled viewer assets used by the extension
  scripts/bundle-viewer.cjs # Regenerates self-contained template bundle
  package.json              # Commands, activation, build/package scripts
```

## WHERE TO LOOK
| Task | Location | Notes |
| --- | --- | --- |
| Editor lifecycle | `src/extension.ts` | `AseviewEditorProvider` registers and refreshes views. |
| File parsing | `node/parse_with_ase_ts.mjs` | Reads all frames, normalizes arrays, fixed atoms, metadata. |
| Copy structure | `node/write_with_ase_ts.mjs` | Supports `xyz`, `extxyz`, `cif`, `vasp`/`POSCAR`. |
| Bundle generation | `scripts/bundle-viewer.cjs` | Inlines media dependencies into `media/molecular_viewer_bundle.html`. |
| User settings | `package.json` | `aseview.defaultViewer`, style, index, format. |

## CONVENTIONS
- The extension does not import Python or ASE; it shells Node workers using
  `process.execPath`.
- Viewer frames must stay schema-compatible with Python `MolecularData` output.
- `media/molecular_viewer.html` and `media/styles.js` are copies of package
  source files; regenerate them rather than hand-editing extension-only drift.
- The webview receives the large self-contained template bundle through
  `postMessage` to avoid embedding a huge string inline.
- Use `retainContextWhenHidden: true` because viewer state is interactive.

## ANTI-PATTERNS
- Do not change worker output fields without checking Python template consumers.
- Do not edit `media/molecular_viewer_bundle.html` directly except for generated
  output review.
- Do not assume single-frame files; `readIndex` can select slices from all
  parsed frames.
- Do not rely on CDN-only behavior for VS Code. Local media/bundle fallbacks are
  part of the extension's offline path.

## COMMANDS
```bash
npm install
npm run compile
npm run build-worker
npm run bundle-viewer
npm run build
npm run package
```

## GOTCHAS
- `extension.ts` contains several validation/coercion helpers near the bottom;
  keep them strict because webview messages are untyped.
- Parse worker timeout is 45 seconds. Large trajectory support should respect
  that path or make timeout behavior explicit.
- Bundle sync is also enforced by `.github/workflows/sync-vscode-bundle.yml`
  and the optional `.github/hooks/pre-commit` hook.
