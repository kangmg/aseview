"""Molecular template camera preset and export bridge contracts."""

from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_PATHS = [
    ROOT / "aseview" / "templates" / "molecular_viewer.html",
    *sorted((ROOT / "aseview" / "themes").glob("*/molecular_viewer.html")),
]


@pytest.fixture(params=TEMPLATE_PATHS, ids=lambda path: str(path.relative_to(ROOT)))
def molecular_template(request):
    return request.param


def _content(path):
    return path.read_text(encoding="utf-8")


def _function_block(content, start, end):
    start_index = content.index(start)
    end_index = content.index(end, start_index)
    return content[start_index:end_index]


def test_existing_ui_exports_and_diagonal_fit_are_present(molecular_template):
    content = _content(molecular_template)

    assert "function saveAsPNG(" in content
    assert "function saveAsGIF(" in content
    assert "fitCameraToMolecule" in content
    assert "new THREE.Vector3(distance, distance, distance)" in content


def test_camera_preset_helpers_and_settings_are_present(molecular_template):
    content = _content(molecular_template)

    required_snippets = [
        "const CARTESIAN_VIEW_DIRECTIONS = Object.freeze({",
        "top: [0, 0, 1]",
        "bottom: [0, 0, -1]",
        "front: [0, 1, 0]",
        "back: [0, -1, 0]",
        "right: [1, 0, 0]",
        "left: [-1, 0, 0]",
        "function isValidCellMatrix(cell)",
        "function resolveCellViewDirection(preset, frame)",
        "function vectorFromEulerXYZ(eulerDegrees)",
        "function resolveViewDirection(viewSpec, frame, strict)",
        "function stableUpVector(direction, requestedUp)",
        "function applyViewSpec(viewSpec, options = {})",
        "function resetView()",
        "viewPreset: null",
        "viewDirection: null",
        "viewEuler: null",
        "viewUp: null",
    ]
    for snippet in required_snippets:
        assert snippet in content, f"{molecular_template} missing {snippet!r}"


def test_gif_export_honors_sizing_options(molecular_template):
    content = _content(molecular_template)
    gif_block = _function_block(content, "function exportGIF(options = {})", "function saveAsPNG()")

    required_snippets = [
        "const originalWidth = canvas.width;",
        "const originalHeight = canvas.height;",
        "const originalBackground = scene.background;",
        "renderer.getClearColor(originalClearColor);",
        "const targetWidth = opts.width || Math.max(1, Math.round(originalWidth * opts.scale));",
        "const targetHeight = opts.height || Math.max(1, Math.round(originalHeight * opts.scale));",
        "renderer.setSize(targetWidth, targetHeight, false);",
        "const transparent = opts.transparent !== undefined ? !!opts.transparent : true;",
        "const bgColor = opts.backgroundColor || settings.backgroundColor || '#111827';",
        "scene.background = originalBackground;",
        "renderer.setClearColor(originalClearColor, originalClearAlpha);",
        "gifWidth: targetWidth",
        "gifHeight: targetHeight",
        "const result = { ok: true, type: 'gif', filename, width: targetWidth, height: targetHeight };",
    ]
    for snippet in required_snippets:
        assert snippet in gif_block, f"{molecular_template} GIF export missing {snippet!r}"


def test_export_command_bridge_contract_is_present(molecular_template):
    content = _content(molecular_template)

    required_snippets = [
        "function buildExportFilename(type, explicitFilename)",
        "function downloadDataUrl(dataUrl, filename)",
        "function exportPNG(options = {})",
        "function exportGIF(options = {})",
        "function saveAsPNG()",
        "function saveAsGIF()",
        "returnDataUrl",
        "backgroundColor",
        "transparent",
        "sampleInterval",
        "case 'setView':",
        "case 'resetView':",
        "case 'exportImage':",
        "postCommandResult(event, msg.requestId,",
        "postCommandError(event, msg.requestId,",
        "type: 'exportResult'",
        "type: 'exportError'",
    ]
    for snippet in required_snippets:
        assert snippet in content, f"{molecular_template} missing {snippet!r}"
