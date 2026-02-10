#!/usr/bin/env python3
"""Generate example viewer HTML files for documentation."""
import os
import sys
import numpy as np

# Add the package to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ase import Atoms
from ase.build import molecule, bulk, fcc111, graphene_nanoribbon, nanotube
from ase.lattice.cubic import FaceCenteredCubic
from aseview import MolecularViewer, OverlayViewer, NormalViewer

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "assets", "viewers")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def create_water_viewer():
    """Create water molecule viewer."""
    water = molecule("H2O")
    viewer = MolecularViewer(water, style="cartoon")
    viewer.save_html(os.path.join(OUTPUT_DIR, "water.html"))
    print("Created water.html")


def create_ethanol_viewer():
    """Create ethanol molecule viewer with neon style."""
    ethanol = molecule("CH3CH2OH")
    viewer = MolecularViewer(ethanol, style="neon", backgroundColor="#000000")
    viewer.save_html(os.path.join(OUTPUT_DIR, "ethanol_neon.html"))
    print("Created ethanol_neon.html")


def create_benzene_viewer():
    """Create benzene molecule viewer with glossy style."""
    benzene = molecule("C6H6")
    viewer = MolecularViewer(benzene, style="glossy")
    viewer.save_html(os.path.join(OUTPUT_DIR, "benzene_glossy.html"))
    print("Created benzene_glossy.html")


def create_trajectory_viewer():
    """Create a simple trajectory animation."""
    # Create a simple vibrating water molecule trajectory
    water = molecule("H2O")
    trajectory = []

    for i in range(20):
        atoms = water.copy()
        # Simple harmonic oscillation
        t = i / 20 * 2 * np.pi
        displacement = 0.1 * np.sin(t)
        positions = atoms.get_positions()
        positions[1, 0] += displacement  # Move H atom
        positions[2, 0] -= displacement  # Move other H atom
        atoms.set_positions(positions)
        # Add fake energy
        atoms.info['energy'] = -10.0 + 0.5 * np.sin(t)
        trajectory.append(atoms)

    viewer = MolecularViewer(trajectory, showEnergyPlot=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "trajectory.html"))
    print("Created trajectory.html")


def create_overlay_viewer():
    """Create overlay comparison of conformers."""
    # Create slightly different conformations
    mol1 = molecule("CH3CH2OH")
    mol2 = mol1.copy()
    mol3 = mol1.copy()

    # Rotate mol2 slightly
    positions2 = mol2.get_positions()
    angle = np.pi / 12  # 15 degrees
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    rotation = np.array([[cos_a, -sin_a, 0], [sin_a, cos_a, 0], [0, 0, 1]])
    mol2.set_positions(positions2 @ rotation.T)

    # Rotate mol3 differently
    positions3 = mol3.get_positions()
    angle = -np.pi / 12
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    rotation = np.array([[cos_a, -sin_a, 0], [sin_a, cos_a, 0], [0, 0, 1]])
    mol3.set_positions(positions3 @ rotation.T)

    viewer = OverlayViewer([mol1, mol2, mol3], colorBy="Molecule", centerMolecules=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "overlay_conformers.html"))
    print("Created overlay_conformers.html")


def create_overlay_colormap_viewer():
    """Create overlay with colormap."""
    # Create optimization-like trajectory
    water = molecule("H2O")
    trajectory = []

    for i in range(10):
        atoms = water.copy()
        # Gradually change geometry
        scale = 1.0 + i * 0.02
        positions = atoms.get_positions()
        positions *= scale
        atoms.set_positions(positions)
        trajectory.append(atoms)

    viewer = OverlayViewer(trajectory, colorBy="Colormap", colormap="viridis", centerMolecules=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "overlay_colormap.html"))
    print("Created overlay_colormap.html")


def create_normal_mode_viewer():
    """Create normal mode viewer with fake vibration data."""
    water = molecule("H2O")

    # Create fake mode vectors (3 atoms × 3 modes for water's real modes)
    # Each mode is displacement per atom [x, y, z]
    n_atoms = len(water)

    # Symmetric stretch
    mode1 = [[0, 0, 0],  # O stays
             [0.7, 0, 0.3],   # H1 moves
             [-0.7, 0, 0.3]]  # H2 moves

    # Asymmetric stretch
    mode2 = [[0, 0, 0.1],
             [0.7, 0, -0.3],
             [0.7, 0, -0.3]]

    # Bending
    mode3 = [[0, 0, -0.2],
             [0, 0.5, 0.4],
             [0, -0.5, 0.4]]

    mode_vectors = [mode1, mode2, mode3]
    frequencies = [1595.0, 3657.0, 3756.0]  # Approximate water frequencies

    viewer = NormalViewer(
        water,
        mode_vectors=mode_vectors,
        frequencies=frequencies,
        showModeVector=True
    )
    viewer.save_html(os.path.join(OUTPUT_DIR, "normal_modes.html"))
    print("Created normal_modes.html")


def create_caffeine_viewer():
    """Create caffeine molecule viewer."""
    # Caffeine structure
    caffeine = Atoms(
        symbols=['C', 'C', 'C', 'C', 'N', 'C', 'N', 'C', 'N', 'C', 'N', 'O', 'O', 'C', 'C', 'C',
                 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
        positions=[
            [0.000, 0.000, 0.000],
            [1.400, 0.000, 0.000],
            [2.100, 1.200, 0.000],
            [1.400, 2.400, 0.000],
            [0.000, 2.400, 0.000],
            [-0.700, 1.200, 0.000],
            [2.100, 3.600, 0.000],
            [3.500, 3.300, 0.000],
            [3.500, 1.800, 0.000],
            [4.700, 0.900, 0.000],
            [-0.700, 3.600, 0.000],
            [-1.900, 1.200, 0.000],
            [4.200, 4.200, 0.000],
            [1.800, 4.900, 0.000],
            [-0.200, 4.900, 0.000],
            [5.100, -0.300, 0.000],
            [2.000, -0.900, 0.000],
            [1.200, 5.100, 0.900],
            [1.200, 5.100, -0.900],
            [2.700, 5.500, 0.000],
            [-0.600, 5.100, 0.900],
            [-0.600, 5.100, -0.900],
            [-1.100, 5.500, 0.000],
            [6.100, -0.100, 0.000],
            [4.900, -0.900, 0.900],
            [4.900, -0.900, -0.900],
        ]
    )
    viewer = MolecularViewer(caffeine, style="metallic")
    viewer.save_html(os.path.join(OUTPUT_DIR, "caffeine.html"))
    print("Created caffeine.html")


def create_silicon_crystal_viewer():
    """Create silicon crystal (diamond structure)."""
    si = bulk('Si', 'diamond', a=5.43, cubic=True)
    # Create 2x2x2 supercell for better visualization
    si = si * (2, 2, 2)
    viewer = MolecularViewer(si, style="glossy", showCell=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "silicon_crystal.html"))
    print("Created silicon_crystal.html")


def create_nacl_crystal_viewer():
    """Create NaCl crystal (rocksalt structure)."""
    nacl = bulk('NaCl', crystalstructure='rocksalt', a=5.64)
    # Create 2x2x2 supercell
    nacl = nacl * (2, 2, 2)
    viewer = MolecularViewer(nacl, style="default", showCell=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "nacl_crystal.html"))
    print("Created nacl_crystal.html")


def create_gold_surface_viewer():
    """Create Au(111) surface slab."""
    # Create 4-layer Au(111) slab with 3x3 surface cell
    au_slab = fcc111('Au', size=(3, 3, 4), vacuum=5.0)
    viewer = MolecularViewer(au_slab, style="metallic", showCell=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "gold_surface.html"))
    print("Created gold_surface.html")


def create_graphene_viewer():
    """Create graphene sheet."""
    # Create graphene nanoribbon (zigzag, width=4, length=6)
    graphene = graphene_nanoribbon(4, 6, type='zigzag', saturated=True, vacuum=5.0)
    viewer = MolecularViewer(graphene, style="cartoon", showCell=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "graphene.html"))
    print("Created graphene.html")


def create_nanotube_viewer():
    """Create carbon nanotube."""
    cnt = nanotube(6, 0, length=4, vacuum=5.0)
    viewer = MolecularViewer(cnt, style="neon", backgroundColor="#000000")
    viewer.save_html(os.path.join(OUTPUT_DIR, "carbon_nanotube.html"))
    print("Created carbon_nanotube.html")


def create_fcc_metal_viewer():
    """Create FCC copper crystal."""
    cu = bulk('Cu', 'fcc', a=3.61, cubic=True)
    cu = cu * (3, 3, 3)
    viewer = MolecularViewer(cu, style="metallic", showCell=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "copper_fcc.html"))
    print("Created copper_fcc.html")


if __name__ == "__main__":
    print("Generating example viewer HTML files...")
    # Molecules
    create_water_viewer()
    create_ethanol_viewer()
    create_benzene_viewer()
    create_caffeine_viewer()
    # Trajectory
    create_trajectory_viewer()
    # Overlay
    create_overlay_viewer()
    create_overlay_colormap_viewer()
    # Normal modes
    create_normal_mode_viewer()
    # Solid-state structures
    create_silicon_crystal_viewer()
    create_nacl_crystal_viewer()
    create_gold_surface_viewer()
    create_graphene_viewer()
    create_nanotube_viewer()
    create_fcc_metal_viewer()
    print("Done!")
