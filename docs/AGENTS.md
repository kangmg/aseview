# DOCS GUIDE

## OVERVIEW
Docs are built with MkDocs and include committed generated viewer HTML examples
under `docs/assets/viewers/` for GitHub Pages.

## STRUCTURE
```
docs/
  index.md
  getting-started/
  examples/
  api/
  assets/viewers/       # Generated HTML viewers
  playground/index.html # Browser playground
  generate_examples.py  # Source of committed viewer HTML examples
```

## WHERE TO LOOK
| Task | Location | Notes |
| --- | --- | --- |
| Navigation/site config | `../mkdocs.yml` | Page order and Material theme config. |
| Example viewer HTML | `generate_examples.py` | Regenerate after viewer/template behavior changes. |
| User-facing CLI docs | `cli-reference.md` | Keep in sync with `aseview/cli.py`. |
| JS module docs | `js-module.md` | Keep in sync with `aseview/static/js/aseview.js`. |
| Theme docs | `theming.md` | Keep in sync with `aseview/themes/` and public theme helpers. |
| API pages | `api/*.md` | Viewer-specific options and behavior. |

## CONVENTIONS
- Generated files in `docs/assets/viewers/` are committed intentionally.
- `docs/generate_examples.py` imports local package code by adding the repo root
  to `sys.path`; run it from the repo root or account for relative paths.
- Docs deploy installs `mkdocs-material`, `pymdown-extensions`, `numpy`, `ase`,
  and `pip install -e .` before generating examples.
- Keep docs command examples aligned with current CLI option names and defaults.

## ANTI-PATTERNS
- Do not hand-edit generated viewer HTML when the source is
  `generate_examples.py` or the viewer templates.
- Do not describe a theme or style unless it exists in package runtime assets.
- Do not let docs examples use APIs that tests or README no longer exercise.
- Default projection is `Perspective`; `Orthographic` remains a supported
  selectable mode in implementation and tests.

## COMMANDS
```bash
python docs/generate_examples.py
mkdocs build
mkdocs serve
```

## GOTCHAS
- Generated viewer HTML is large because dependencies are inlined.
- Pages deployment rebuilds examples on every docs/main package change, so docs
  failures can come from runtime viewer regressions.
