"""Focused tests for Python-side viewer camera settings."""

import json

import pytest
from ase.build import molecule

from aseview.wrapper import MolecularViewer, NormalViewer, OverlayViewer


CAMERA_DEFAULTS = {
    "viewPreset": None,
    "viewDirection": None,
    "viewEuler": None,
    "viewUp": None,
    "viewFit": 1.0,
}


@pytest.fixture
def h2o():
    return molecule("H2O")


@pytest.mark.parametrize(
    "viewer_cls",
    [MolecularViewer, OverlayViewer, NormalViewer],
)
def test_camera_settings_defaults_are_serializable(viewer_cls, h2o):
    viewer = viewer_cls(h2o)

    for key, expected in CAMERA_DEFAULTS.items():
        assert viewer.settings[key] == expected

    assert viewer.settings["viewMode"] == "Orthographic"
    json.dumps(viewer.settings)


@pytest.mark.parametrize(
    "viewer_cls",
    [MolecularViewer, OverlayViewer, NormalViewer],
)
def test_custom_camera_settings_propagate_to_html(viewer_cls, h2o):
    viewer = viewer_cls(
        h2o,
        viewPreset="top-c",
        viewDirection=[0, 0, 1],
        viewUp=[0, 1, 0],
        viewFit=1.25,
    )

    html = viewer.get_html()

    assert '"viewPreset": "top-c"' in html
    assert '"viewDirection": [0, 0, 1]' in html
    assert '"viewUp": [0, 1, 0]' in html
    assert '"viewFit": 1.25' in html
