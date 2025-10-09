#!/usr/bin/env python3
"""
Example usage of aseview molecular viewer.
"""

import numpy as np
from ase import Atoms
from ase.build import molecule

import aseview


def example_molecular_viewer():
    """Example of using the MolecularViewer with ASE Atoms."""
    # Create a water molecule
    atoms = molecule('H2O')
    print(f"Created water molecule with {len(atoms)} atoms")
    
    # Create viewer
    viewer = aseview.MolecularViewer(atoms)
    
    # Save to HTML file
    viewer.save_html('molecular_viewer_example.html')
    print("Saved molecular viewer example to 'molecular_viewer_example.html'")


def example_normal_viewer():
    """Example of using the NormalViewer."""
    # Create a simple trajectory-like data
    atoms1 = molecule('H2O')
    atoms2 = atoms1.copy()
    atoms2.positions += np.array([0.1, 0.0, 0.0])  # Small displacement
    
    # Create viewer with multiple configurations
    viewer = aseview.NormalViewer([atoms1, atoms2])
    
    # Save to HTML file
    viewer.save_html('normal_viewer_example.html')
    print("Saved normal viewer example to 'normal_viewer_example.html'")


def example_overlay_viewer():
    """Example of using the OverlayViewer."""
    # Create different molecules
    water = molecule('H2O')
    ammonia = molecule('NH3')
    
    # Create viewer with multiple molecules
    viewer = aseview.OverlayViewer([water, ammonia])
    
    # Save to HTML file
    viewer.save_html('overlay_viewer_example.html')
    print("Saved overlay viewer example to 'overlay_viewer_example.html'")


def example_from_json():
    """Example of using the viewers with JSON data."""
    # Load data from the existing data.json file
    viewer = aseview.MolecularViewer('data.json')
    
    # Save to HTML file
    viewer.save_html('molecular_viewer_from_json.html')
    print("Saved molecular viewer from JSON example to 'molecular_viewer_from_json.html'")


if __name__ == "__main__":
    print("Creating aseview examples...")
    
    example_molecular_viewer()
    example_normal_viewer()
    example_overlay_viewer()
    example_from_json()
    
    print("Done! Open the HTML files in your browser to view the examples.")