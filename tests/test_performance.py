"""Test performance-related settings and code paths."""

import pytest
from ase import Atoms
from ase.build import molecule

from aseview.wrapper import MolecularViewer, NormalViewer, OverlayViewer


def _make_large_structure(n_atoms=500):
    """Create a large structure for performance testing."""
    import numpy as np

    rng = np.random.default_rng(42)
    symbols = ["C"] * n_atoms
    positions = rng.uniform(-20, 20, size=(n_atoms, 3))
    return Atoms(symbols=symbols, positions=positions)


class TestPerformanceMode:
    def test_auto_mode_default(self):
        viewer = MolecularViewer(molecule("H2O"))
        assert viewer.settings["performanceMode"] == "auto"

    def test_high_mode(self):
        viewer = MolecularViewer(molecule("H2O"), performanceMode="high")
        assert viewer.settings["performanceMode"] == "high"

    def test_normal_mode(self):
        viewer = MolecularViewer(molecule("H2O"), performanceMode="normal")
        assert viewer.settings["performanceMode"] == "normal"

    def test_performance_mode_in_html(self):
        viewer = MolecularViewer(molecule("H2O"), performanceMode="high")
        html = viewer.get_html()
        assert '"performanceMode"' in html or "'performanceMode'" in html

    def test_normal_viewer_has_performance_mode(self):
        viewer = NormalViewer(molecule("H2O"))
        assert "performanceMode" in viewer.settings

    def test_overlay_viewer_has_performance_mode(self):
        viewer = OverlayViewer(molecule("H2O"))
        assert "performanceMode" in viewer.settings


class TestLargeStructureHTML:
    """Ensure HTML generation works for larger structures without errors."""

    def test_large_structure_html_generation(self):
        atoms = _make_large_structure(200)
        viewer = MolecularViewer(atoms)
        html = viewer.get_html()
        assert isinstance(html, str)
        assert len(html) > 1000

    def test_large_structure_with_high_perf(self):
        atoms = _make_large_structure(200)
        viewer = MolecularViewer(atoms, performanceMode="high")
        html = viewer.get_html()
        assert "getLODSegments" in html

    def test_html_contains_lod_functions(self):
        atoms = _make_large_structure(50)
        viewer = MolecularViewer(atoms)
        html = viewer.get_html()
        assert "getLODSegments" in html
        assert "getCachedSphereGeometry" in html
        assert "getCachedMaterial" in html
        assert "findBondsSpatialGrid" in html
        assert "clearPerformanceCache" in html
