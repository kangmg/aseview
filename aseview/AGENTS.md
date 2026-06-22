# ASEVIEW PACKAGE GUIDE

## OVERVIEW
This directory is the runtime Python package plus the shipped HTML/JS assets
that render structures in notebooks, saved HTML, the browser module, and docs.

## STRUCTURE
```
aseview/
  wrapper.py            # Core data conversion and viewer classes
  cli.py                # Typer CLI and local HTML server
  hessian_parsers.py    # ORCA/VASP normal-mode parsers
  jupyter.py            # ipywidgets wrapper
  static/js/            # Browser API and shared render helpers
  templates/            # Backward-compatible viewer templates
  themes/<name>/        # Complete themed template sets
```

## WHERE TO LOOK
| Task | Location | Notes |
| --- | --- | --- |
| Add viewer option | `wrapper.py` | Add Python setting, then confirm template consumes it. |
| Add atom/frame metadata | `MolecularData.from_atoms()` | Keep VS Code `parse_with_ase_ts.mjs` schema parity in mind. |
| Theme resolution | `_resolve_template()` | Searches `themes/<theme>/` before `templates/`. |
| Shared render style | `static/js/styles.js` | Atom/bond factories and performance utilities. |
| CLI option | `cli.py` | Typer option plus viewer kwargs wiring. |
| Normal-mode import | `hessian_parsers.py` and `NormalViewer.from_*()` | Preserve imaginary-mode handling. |

## CONVENTIONS
- Public API names are re-exported from `__init__.py`; add new user-facing
  classes/helpers there deliberately.
- Viewer HTML generation inlines local JS dependencies so notebook output and
  saved HTML can work without a live CDN.
- Packaging asset changes need both `pyproject.toml` package-data and
  `MANIFEST.in` considered; they do not currently list identical asset sets.
- Settings are plain dictionaries merged with `**kwargs`; template-side names
  should match the Python keys exactly.
- Theme directories are full template sets, not partial overrides.
- `static/js/aseview.js` loads templates into iframes and talks to them through
  postMessage; it is not the renderer implementation.

## ANTI-PATTERNS
- Do not introduce a template placeholder in one viewer path without checking
  Python, browser JS, docs examples, and VS Code use.
- Do not leave only `templates/` updated when themed files need identical
  behavior.
- Do not assume all input is a single `Atoms`; `BaseViewer._process_data()`
  accepts paths, globs, lists, dicts, and iterable trajectories.
- Do not remove vendor files from package data; they are used for inlining.

## CHECKS
```bash
python -m compileall aseview/ tests/ -q
pytest tests/test_imports.py tests/test_viewers.py tests/test_molecular_data.py -v
pytest tests/test_js_syntax.py -v
node --check aseview/static/js/styles.js
node --check aseview/static/js/aseview.js
```

## GOTCHAS
- `MolecularViewer`, `NormalViewer`, `OverlayViewer`, `FragSelector`, and
  `LiteViewer` each load different templates but share the same data contract.
- `NormalViewer.from_vasp()` may need to expand selective-dynamics modes back
  onto the full structure.
- `LiteViewer` has alias parameters like `styles`, `hide_hs`, and `centering`
  that tests cover.
