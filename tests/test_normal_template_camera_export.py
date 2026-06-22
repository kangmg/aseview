from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
BASE_NORMAL_TEMPLATE = ROOT / "aseview" / "templates" / "normal_viewer.html"
THEMES_DIR = ROOT / "aseview" / "themes"


def _normal_template_paths():
    paths = [BASE_NORMAL_TEMPLATE]
    paths.extend(sorted(THEMES_DIR.glob("*/normal_viewer.html")))
    return paths


def _function_block(content, start, end):
    start_index = content.index(start)
    end_index = content.index(end, start_index)
    return content[start_index:end_index]


@pytest.fixture(params=_normal_template_paths(), ids=lambda path: str(path.relative_to(ROOT)))
def normal_template(request):
    return request.param


def test_existing_normal_export_fit_and_projection_baseline(normal_template):
    content = normal_template.read_text(encoding="utf-8")

    assert "function saveAsPNG()" in content
    assert "function saveAsGIF()" in content
    assert "function fitCameraToMolecule(group)" in content
    assert "new THREE.Vector3(distance, distance, distance)" in content
    assert "viewMode: 'Perspective'" in content


def test_normal_templates_expose_camera_and_export_command_contract(normal_template):
    content = normal_template.read_text(encoding="utf-8")

    required_snippets = [
        "viewPreset: null",
        "viewDirection: null",
        "viewEuler: null",
        "viewUp: null",
        "viewFit: 1",
        "function getNormalModeBaseFrame()",
        "function isValidCellMatrix(cell)",
        "function resolveViewSpec(viewSpec, options = {})",
        "function applyCameraView(viewSpec, options = {})",
        "function applyInitialViewSettings()",
        "function resetView()",
        "function exportNormalPNG(options = {})",
        "function exportNormalGIF(options = {})",
        "function saveAsPNG(options)",
        "function saveAsGIF(options)",
        "case 'setView':",
        "case 'resetView':",
        "case 'exportImage':",
        "type: 'exportResult'",
        "code: 'invalid_view_direction'",
        "code: 'missing_vibration_data'",
        "sampleInterval",
    ]

    for snippet in required_snippets:
        assert snippet in content, f"{normal_template} missing {snippet!r}"

    for preset in ["top-c", "bottom-c", "side-a", "side-b", "front", "back", "right", "left"]:
        assert preset in content


def test_normal_gif_export_honors_sizing_options(normal_template):
    content = normal_template.read_text(encoding="utf-8")
    gif_block = _function_block(content, "function exportGIF(options = {})", "function saveAsPNG()")

    required_snippets = [
        "const originalWidth = canvas.width;",
        "const originalHeight = canvas.height;",
        "const targetWidth = opts.width || Math.max(1, Math.round(originalWidth * opts.scale));",
        "const targetHeight = opts.height || Math.max(1, Math.round(originalHeight * opts.scale));",
        "renderer.setSize(targetWidth, targetHeight, false);",
        "gifWidth: targetWidth",
        "gifHeight: targetHeight",
        "const result = { ok: true, type: 'gif', filename, width: targetWidth, height: targetHeight };",
    ]

    for snippet in required_snippets:
        assert snippet in gif_block, f"{normal_template} GIF export missing {snippet!r}"
