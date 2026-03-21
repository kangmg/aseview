"""Test that all package modules import without errors."""

import pytest


def test_import_package():
    """Package-level import should succeed."""
    import aseview
    assert hasattr(aseview, "MolecularViewer")
    assert hasattr(aseview, "NormalViewer")
    assert hasattr(aseview, "OverlayViewer")
    assert hasattr(aseview, "FragSelector")
    assert hasattr(aseview, "MolecularData")


def test_import_wrapper():
    """Wrapper module should import cleanly."""
    from aseview.wrapper import (
        BaseViewer,
        MolecularViewer,
        NormalViewer,
        OverlayViewer,
        FragSelector,
        MolecularData,
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
