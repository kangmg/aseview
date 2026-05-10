"""Test that all package modules import without errors."""

import pytest


def test_import_package():
    """Package-level import should succeed."""
    import aseview
    assert hasattr(aseview, "MolecularViewer")
    assert hasattr(aseview, "NormalViewer")
    assert hasattr(aseview, "OverlayViewer")
    assert hasattr(aseview, "FragSelector")
    assert hasattr(aseview, "LiteViewer")
    assert hasattr(aseview, "view")
    assert hasattr(aseview, "MolecularData")
    assert hasattr(aseview, "set_theme")
    assert hasattr(aseview, "get_theme")
    assert hasattr(aseview, "list_themes")


def test_theme_api():
    """Theme management functions should work correctly."""
    import aseview
    original = aseview.get_theme()
    assert isinstance(original, str)

    themes = aseview.list_themes()
    assert isinstance(themes, list)
    assert "dark" in themes

    aseview.set_theme("dark")
    assert aseview.get_theme() == "dark"

    aseview.set_theme(original)  # restore


def test_import_wrapper():
    """Wrapper module should import cleanly."""
    from aseview.wrapper import (
        BaseViewer,
        MolecularViewer,
        NormalViewer,
        OverlayViewer,
        FragSelector,
        LiteViewer,
        view,
        MolecularData,
        set_theme,
        get_theme,
        list_themes,
        _resolve_template,
    )


def test_import_utils():
    """Utils module should import cleanly."""
    from aseview import utils


def test_import_cli():
    """CLI module should import cleanly."""
    from aseview import cli


def test_import_jupyter():
    """Jupyter module should import cleanly."""
    from aseview import jupyter


def test_import_hessian_parsers():
    """Hessian parsers module should import cleanly."""
    from aseview import hessian_parsers
