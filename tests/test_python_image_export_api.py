import base64

import pytest
from ase import Atoms

import aseview.export as export_module
from aseview import ImageExportError, MolecularViewer, OverlayViewer


PNG_DATA_URL = "data:image/png;base64," + base64.b64encode(b"png-bytes").decode("ascii")
GIF_DATA_URL = "data:image/gif;base64," + base64.b64encode(b"gif-bytes").decode("ascii")


def _single_atom_viewer():
    return MolecularViewer(Atoms("H", positions=[[0, 0, 0]]), viewPreset="top-c")


def test_save_png_writes_file_and_forwards_render_options(monkeypatch, tmp_path):
    captured = {}

    def fake_render(html, *, export_format, options, timeout_ms, browser_executable_path):
        captured.update(
            html=html,
            export_format=export_format,
            options=options,
            timeout_ms=timeout_ms,
            browser_executable_path=browser_executable_path,
        )
        return {
            "ok": True,
            "type": "png",
            "filename": options["filename"],
            "width": options["width"],
            "height": options["height"],
            "dataUrl": PNG_DATA_URL,
        }

    monkeypatch.setattr(export_module, "_render_export_data_url", fake_render)

    output = tmp_path / "molecule.png"
    result = _single_atom_viewer().save_png(
        str(output),
        width=1200,
        height=900,
        scale=2,
        transparent=False,
        background_color="#ffffff",
        timeout_ms=1234,
        browser_executable_path="/usr/bin/chromium",
    )

    assert output.read_bytes() == b"png-bytes"
    assert result["path"] == str(output)
    assert result["bytes"] == len(b"png-bytes")
    assert captured["export_format"] == "png"
    assert captured["timeout_ms"] == 1234
    assert captured["browser_executable_path"] == "/usr/bin/chromium"
    assert "viewPreset" in captured["html"]
    assert captured["options"] == {
        "download": False,
        "returnDataUrl": True,
        "filename": "molecule.png",
        "scale": 2.0,
        "transparent": False,
        "width": 1200,
        "height": 900,
        "backgroundColor": "#ffffff",
    }


def test_save_png_quality_maps_to_resolution_scale(monkeypatch, tmp_path):
    captured = {}

    def fake_render(html, *, export_format, options, timeout_ms, browser_executable_path):
        captured["options"] = options
        return {
            "ok": True,
            "type": "png",
            "filename": options["filename"],
            "dataUrl": PNG_DATA_URL,
        }

    monkeypatch.setattr(export_module, "_render_export_data_url", fake_render)

    _single_atom_viewer().save_png(str(tmp_path / "high.png"), quality="high")

    assert captured["options"]["scale"] == 3.0


def test_save_gif_writes_file_and_maps_quality(monkeypatch, tmp_path):
    captured = {}

    def fake_render(html, *, export_format, options, timeout_ms, browser_executable_path):
        captured["export_format"] = export_format
        captured["options"] = options
        return {
            "ok": True,
            "type": "gif",
            "filename": options["filename"],
            "dataUrl": GIF_DATA_URL,
        }

    monkeypatch.setattr(export_module, "_render_export_data_url", fake_render)

    output = tmp_path / "mode.gif"
    result = _single_atom_viewer().save_gif(
        str(output),
        frames=12,
        fps=24,
        delay=40,
        quality="high",
        transparent=False,
        background_color="#111827",
    )

    assert output.read_bytes() == b"gif-bytes"
    assert result["path"] == str(output)
    assert captured["export_format"] == "gif"
    assert captured["options"]["frames"] == 12
    assert captured["options"]["fps"] == 24.0
    assert captured["options"]["delay"] == 40.0
    assert captured["options"]["sampleInterval"] == 5
    assert captured["options"]["transparent"] is False
    assert captured["options"]["backgroundColor"] == "#111827"


def test_save_gif_sample_interval_overrides_quality(monkeypatch, tmp_path):
    captured = {}

    def fake_render(html, *, export_format, options, timeout_ms, browser_executable_path):
        captured["options"] = options
        return {
            "ok": True,
            "type": "gif",
            "filename": options["filename"],
            "dataUrl": GIF_DATA_URL,
        }

    monkeypatch.setattr(export_module, "_render_export_data_url", fake_render)

    _single_atom_viewer().save_gif(
        str(tmp_path / "mode.gif"),
        sample_interval=7,
        quality="high",
    )

    assert captured["options"]["sampleInterval"] == 7


def test_save_png_and_save_gif_validate_quality(tmp_path):
    viewer = _single_atom_viewer()

    with pytest.raises(ValueError, match="PNG quality"):
        viewer.save_png(str(tmp_path / "bad.png"), quality="ultra")

    with pytest.raises(ValueError, match="GIF quality"):
        viewer.save_gif(str(tmp_path / "bad.gif"), quality="ultra")

    with pytest.raises(ValueError, match="sample_interval"):
        viewer.save_gif(str(tmp_path / "bad.gif"), sample_interval=0)


def test_decode_rejects_wrong_export_data_url():
    with pytest.raises(ImageExportError, match="invalid PNG data"):
        export_module._decode_image_data_url(GIF_DATA_URL, "png")


def test_overlay_viewer_exposes_png_but_rejects_gif_in_template(monkeypatch, tmp_path):
    captured = {}

    def fake_render(html, *, export_format, options, timeout_ms, browser_executable_path):
        captured["html"] = html
        captured["export_format"] = export_format
        return {
            "ok": True,
            "type": "png",
            "filename": options["filename"],
            "dataUrl": PNG_DATA_URL,
        }

    monkeypatch.setattr(export_module, "_render_export_data_url", fake_render)
    atoms = Atoms("H", positions=[[0, 0, 0]])
    viewer = OverlayViewer([atoms, atoms.copy()])

    viewer.save_png(str(tmp_path / "overlay.png"))

    assert captured["export_format"] == "png"
    assert "Overlay GIF export is not supported" in captured["html"]
