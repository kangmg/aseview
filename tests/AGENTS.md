# TESTS GUIDE

## OVERVIEW
Tests cover import surfaces, viewer settings and HTML generation, frame-schema
conversion, performance helpers, and template/JS syntax invariants.

## STRUCTURE
```
tests/
  test_imports.py        # Package and optional module imports
  test_molecular_data.py # ASE Atoms <-> viewer-frame schema
  test_viewers.py        # Viewer constructors, settings, themes, helpers
  test_performance.py    # Large-structure/performance-mode expectations
  test_js_syntax.py      # JS files and all template script blocks
```

## WHERE TO LOOK
| Task | Location | Notes |
| --- | --- | --- |
| New public export | `test_imports.py` | Add smoke coverage for import behavior. |
| New frame field | `test_molecular_data.py` | Assert Python conversion and roundtrip behavior. |
| New viewer setting | `test_viewers.py` | Check default and custom settings. |
| Template behavior | `test_js_syntax.py` | Add snippets that must stay wired across themes. |
| Large rendering path | `test_performance.py` | Guard LOD/performance utility presence. |

## CONVENTIONS
- Tests intentionally inspect generated HTML strings because browser automation
  is not part of the current suite.
- `test_js_syntax.py` runs `node --check` when Node exists and falls back to a
  bracket-balance heuristic for inline template scripts.
- Theme tests iterate across `aseview/themes/*` plus fallback templates.
- Keep assertions on concrete snippets when a template function or UI hook is a
  cross-surface contract.

## ANTI-PATTERNS
- Do not add tests that require opening a browser or VS Code unless the project
  gains an explicit E2E harness.
- Do not skip themed templates when changing base template behavior.
- Do not test only `MolecularViewer` when the schema or setting affects
  overlay, normal mode, fragment selector, or lite viewer paths.

## COMMANDS
```bash
pytest tests/ -v --tb=short
pytest tests/test_js_syntax.py -v
python -m compileall aseview/ tests/ -q
node --check aseview/static/js/styles.js
node --check aseview/static/js/aseview.js
```

## GOTCHAS
- Node-dependent syntax checks are skipped only when Node is unavailable.
- Template tests include behavior alignment checks such as perspective defaults,
  frag-selector style options, and fixed-constraint highlighting wiring.
