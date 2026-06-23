from pathlib import Path
import re

import pytest


ROOT = Path(__file__).resolve().parents[1]
OVERLAY_TEMPLATES = [
    ROOT / "aseview" / "templates" / "overlay_viewer.html",
    *sorted((ROOT / "aseview" / "themes").glob("*/overlay_viewer.html")),
]


@pytest.mark.parametrize("template_path", OVERLAY_TEMPLATES, ids=lambda p: str(p.relative_to(ROOT)))
def test_overlay_template_camera_export_contract(template_path):
    content = template_path.read_text(encoding="utf-8")

    baseline_snippets = [
        "viewMode: 'Orthographic'",
        "function saveAsPNG(",
        "function fitCameraToMolecule(group)",
    ]
    for snippet in baseline_snippets:
        assert snippet in content, f"{template_path} missing legacy marker {snippet!r}"

    required_snippets = [
        "viewPreset",
        "viewDirection",
        "viewEuler",
        "viewUp",
        "viewFit",
        "function resolveViewDirection(",
        "function applyViewSpec(",
        "function resetView(",
        "function exportPNG(",
        "function saveAsGIF(",
        "case 'setView':",
        "case 'resetView':",
        "case 'exportImage':",
        "unsupported_export",
        "invalid_view_direction",
        "top-c",
        "bottom-c",
        "side-a",
        "side-b",
        "returnDataUrl",
        "backgroundColor",
        "transparent",
        "event.source.postMessage",
        "type: 'png'",
    ]
    for snippet in required_snippets:
        assert snippet in content, f"{template_path} missing overlay camera/export contract {snippet!r}"

    active_view_declarations = re.findall(r"\b(const|let|var)\s+activeViewSpec\b", content)
    assert active_view_declarations == ["let"], (
        f"{template_path} must declare activeViewSpec exactly once; "
        f"found {len(active_view_declarations)} declarations"
    )


@pytest.mark.parametrize("template_path", OVERLAY_TEMPLATES, ids=lambda p: str(p.relative_to(ROOT)))
def test_overlay_template_rejects_gif_export(template_path):
    content = template_path.read_text(encoding="utf-8")

    assert "type: 'gif'" in content
    assert "unsupported_export" in content
    assert "Overlay GIF export is not supported" in content
