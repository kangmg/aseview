"""
aseview - A molecular viewer for ASE (Atomic Simulation Environment) data
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("aseview")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = "aseview"

from .wrapper import MolecularViewer, NormalViewer, OverlayViewer, MolecularData

__all__ = ["MolecularViewer", "NormalViewer", "OverlayViewer", "MolecularData"]
