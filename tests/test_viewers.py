"""Test viewer instantiation and HTML generation."""

import os
import pytest
from ase import Atoms
from ase.build import molecule, bulk

from aseview.wrapper import (
    MolecularViewer, NormalViewer, OverlayViewer, FragSelector,
    set_theme, get_theme, list_themes, _resolve_template,
)


@pytest.fixture
def h2o():
    return molecule("H2O")


@pytest.fixture
def ethanol():
    return molecule("CH3CH2OH")


@pytest.fixture
def si_bulk():
    return bulk("Si", "diamond", a=5.43)


# ── MolecularViewer ──────────────────────────────────────────


class TestMolecularViewer:
    def test_init_from_atoms(self, h2o):
        viewer = MolecularViewer(h2o)
        assert viewer.data is not None

    def test_init_from_dict(self, h2o):
        from aseview.wrapper import MolecularData
        data = MolecularData.from_atoms(h2o)
        viewer = MolecularViewer(data)
        assert viewer.data is not None

    def test_init_from_list(self, h2o, ethanol):
        viewer = MolecularViewer([h2o, ethanol])
        assert isinstance(viewer.data, list)
        assert len(viewer.data) == 2

    def test_default_settings(self, h2o):
        viewer = MolecularViewer(h2o)
        assert viewer.settings["style"] == "cartoon"
        assert viewer.settings["performanceMode"] == "auto"
        assert viewer.settings["showBond"] is True

    def test_custom_settings(self, h2o):
        viewer = MolecularViewer(h2o, style="glossy", performanceMode="high")
        assert viewer.settings["style"] == "glossy"
        assert viewer.settings["performanceMode"] == "high"

    def test_generate_html(self, h2o):
        viewer = MolecularViewer(h2o)
        html = viewer.get_html()
        assert isinstance(html, str)
        assert len(html) > 1000
        assert "THREE" in html

    def test_html_contains_molecular_data(self, h2o):
        viewer = MolecularViewer(h2o)
        html = viewer.get_html()
        # Should contain atom symbols embedded as JSON
        assert "H" in html
        assert "O" in html

    def test_html_contains_performance_utils(self, h2o):
        viewer = MolecularViewer(h2o)
        html = viewer.get_html()
        assert "getLODSegments" in html
        assert "findBondsSpatialGrid" in html
        assert "clearPerformanceCache" in html

    def test_periodic_structure(self, si_bulk):
        viewer = MolecularViewer(si_bulk)
        html = viewer.get_html()
        assert isinstance(html, str)
        assert len(html) > 1000

    def test_save_html(self, h2o, tmp_path):
        viewer = MolecularViewer(h2o)
        filepath = tmp_path / "test.html"
        viewer.save_html(str(filepath))
        assert filepath.exists()
        content = filepath.read_text()
        assert "THREE" in content


# ── NormalViewer ─────────────────────────────────────────────


class TestNormalViewer:
    def test_init_without_vibrations(self, h2o):
        viewer = NormalViewer(h2o)
        assert viewer.mode_vectors is None
        assert viewer.frequencies is None

    def test_default_settings(self, h2o):
        viewer = NormalViewer(h2o)
        assert viewer.settings["performanceMode"] == "auto"
        assert viewer.settings["style"] == "cartoon"

    def test_generate_html(self, h2o):
        viewer = NormalViewer(h2o)
        html = viewer.get_html()
        assert isinstance(html, str)
        assert len(html) > 1000


# ── OverlayViewer ────────────────────────────────────────────


class TestOverlayViewer:
    def test_init_single(self, h2o):
        viewer = OverlayViewer(h2o)
        assert isinstance(viewer.data, list)

    def test_init_multiple(self, h2o, ethanol):
        viewer = OverlayViewer([h2o, ethanol])
        assert len(viewer.data) == 2

    def test_default_settings(self, h2o):
        viewer = OverlayViewer(h2o)
        assert viewer.settings["performanceMode"] == "auto"

    def test_generate_html(self, h2o):
        viewer = OverlayViewer(h2o)
        html = viewer.get_html()
        assert isinstance(html, str)
        assert "findBondsSpatialGrid" in html

    def test_index_list_same_mode(self, h2o):
        viewer = OverlayViewer([h2o, h2o], index_list=[0, 1])
        for mol in viewer.data:
            assert len(mol["symbols"]) == 2

    def test_index_list_per_structure(self, h2o):
        viewer = OverlayViewer([h2o, h2o], index_list=[[0, 1], [1, 2]])
        assert len(viewer.data[0]["symbols"]) == 2
        assert len(viewer.data[1]["symbols"]) == 2


# ── FragSelector ─────────────────────────────────────────────


class TestFragSelector:
    def test_init(self, h2o):
        viewer = FragSelector(h2o)
        assert viewer.data is not None

    def test_generate_html(self, h2o):
        viewer = FragSelector(h2o)
        html = viewer.get_html()
        assert isinstance(html, str)
        assert len(html) > 1000


# ── Theme System ─────────────────────────────────────────────


class TestThemeSystem:
    def test_list_themes(self):
        themes = list_themes()
        assert isinstance(themes, list)
        assert "dark" in themes

    def test_get_set_theme(self):
        original = get_theme()
        set_theme("dark")
        assert get_theme() == "dark"
        set_theme(original)

    def test_resolve_template_dark(self):
        path = _resolve_template("molecular_viewer.html", "dark")
        assert os.path.exists(path)
        assert "dark" in path

    def test_resolve_template_fallback(self):
        path = _resolve_template("molecular_viewer.html", "nonexistent_theme_xyz")
        assert os.path.exists(path)

    def test_resolve_template_all_themes(self):
        for theme in list_themes():
            for name in ["molecular_viewer.html", "normal_viewer.html",
                         "overlay_viewer.html", "frag_selector.html"]:
                path = _resolve_template(name, theme)
                assert os.path.exists(path), f"Missing {theme}/{name}"

    def test_viewer_theme_param(self, h2o):
        for theme in list_themes():
            v = MolecularViewer(h2o, theme=theme)
            html = v.get_html()
            assert len(html) > 1000

    def test_overlay_viewer_theme_param(self, h2o):
        for theme in list_themes():
            v = OverlayViewer(h2o, theme=theme)
            html = v.get_html()
            assert len(html) > 1000

    def test_frag_selector_theme_param(self, h2o):
        for theme in list_themes():
            v = FragSelector(h2o, theme=theme)
            html = v.get_html()
            assert len(html) > 1000

    def test_global_theme_affects_viewers(self, h2o):
        original = get_theme()
        try:
            set_theme("dark")
            v = MolecularViewer(h2o)
            html = v.get_html()
            assert len(html) > 1000
        finally:
            set_theme(original)
