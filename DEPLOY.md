# Deployment Guide for Agents

This repository has two release surfaces:

- Python package: `aseview` on PyPI, version in `pyproject.toml`
- VS Code extension: `aseview-vscode`, version in `vscode-extension/package.json`

Both publish workflows react to `v*` tags. Treat a tag like `v0.0.10` as a shared release unless there is a deliberate reason not to.

Current workflow behavior:

- PyPI publish happens on a `v*` tag or GitHub Release event. Do not rely on `workflow_dispatch target=pypi`; the current publish job does not publish PyPI from that input.
- VSIX asset creation happens on a `v*` tag, GitHub Release event, or `workflow_dispatch` with `release_tag`.
- Visual Studio Marketplace/Open VSX publish only happens from `workflow_dispatch`.

## Golden Rule

Do not create or push a `v*` tag until both package versions and local build checks are finished.

For a normal release, keep these versions identical:

- `pyproject.toml` `[project].version`
- `aseview/static/js/aseview.js` exported/viewer version, if present
- `vscode-extension/package.json` `version`
- `vscode-extension/package-lock.json` root version

If the shared tag is `vX.Y.Z`, all of the above should normally be `X.Y.Z`.

## Recommended Release Flow

1. Pull and inspect state.

   ```bash
   git pull
   git status --short
   ```

2. Choose the next version.

   Never reuse a version that already exists on PyPI, VS Code marketplaces, or GitHub tags. Do not move a tag after PyPI has published from it.

3. Update Python package version.

   Edit:

   - `pyproject.toml`
   - `aseview/static/js/aseview.js`, if its version constant exists
   - `uv.lock`, by running a `uv` command after editing

4. Update VS Code extension version.

   ```bash
   cd vscode-extension
   npm version X.Y.Z --no-git-tag-version
   cd ..
   ```

5. Verify locally.

   ```bash
   uv run --extra dev pytest
   uv run --with build python -m build --sdist --wheel --outdir /private/tmp/aseview-dist-check
   cd vscode-extension
   npm run build
   npm run package -- --out /private/tmp/aseview-vscode-X.Y.Z.vsix
   cd ..
   ```

6. Commit first, then tag.

   ```bash
   git status --short
   git add pyproject.toml uv.lock aseview/static/js/aseview.js vscode-extension/package.json vscode-extension/package-lock.json
   git commit -m "Release aseview X.Y.Z"
   git tag vX.Y.Z
   git push origin main
   git push origin vX.Y.Z
   ```

7. Verify GitHub Actions.

   The pushed tag should trigger:

   - `Publish to PyPI`: builds and publishes `aseview` to PyPI via trusted publisher/OIDC
   - `Publish VS Code Extension`: builds `aseview-vscode-X.Y.Z.vsix` and `aseview-latest.vsix`, then attaches them to the GitHub Release for `vX.Y.Z`

8. Publish VS Code marketplaces only when requested.

   The tag flow creates the VSIX assets, but Marketplace/Open VSX publishing is manual. Run the `Publish VS Code Extension` workflow with:

   - `target`: `marketplace`, `openvsx`, or `both`
   - `release_tag`: `vX.Y.Z` if the GitHub Release asset also needs to be created or refreshed

   Required repository secrets:

   - `VSCE_PAT` for Visual Studio Marketplace
   - `OVSX_PAT` for Open VSX

## Avoiding the Known Failure

Failure:

```text
VSIX artifact not found at dist/aseview-vscode-X.Y.Z.vsix
```

Cause:

- The tag version and VS Code extension package version were different, or
- The workflow expected a filename derived from the tag instead of the actual packaged extension version.

Prevention:

- Bump `vscode-extension/package.json` before tagging.
- Run `npm run package -- --out /private/tmp/aseview-vscode-X.Y.Z.vsix`.
- Confirm the generated VSIX filename uses the intended version.
- The workflow now uses the actual packaged VSIX filename, but agents should still keep the extension version aligned with the tag to avoid surprising release assets.

If the tag workflow already failed:

- Do not move or delete the tag after PyPI has published.
- Fix the workflow/package version on `main`.
- Run `Publish VS Code Extension` manually with `release_tag: vX.Y.Z`.

If a manual run creates the `aseview-vsix` Actions artifact but the GitHub Release is still missing:

- Check the `attach-release-asset` job. If it is `skipped`, the manual run did not receive a non-empty `release_tag`.
- Start a new `Run workflow` execution from `main`; do not use "Re-run jobs".
- Set `target: none` and `release_tag: vX.Y.Z`.

## Trusted Publisher Notes

PyPI publishing uses GitHub trusted publisher/OIDC, not an npm token. Deleting an npm token does not affect PyPI trusted publishing.

Trusted publisher expiration emails mean PyPI is waiting for a release through that configured GitHub workflow/environment. A successful PyPI publish through the tag workflow satisfies it.
