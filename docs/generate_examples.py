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
from aseview import MolecularViewer, OverlayViewer, NormalViewer, FragSelector

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
    mol1 = molecule("CH3CH2OH")

    # Center mol1 at origin so all conformers share the same COM = [0,0,0]
    pos = mol1.get_positions()
    mol1.set_positions(pos - pos.mean(axis=0))

    mol2 = mol1.copy()
    mol3 = mol1.copy()

    mol1.info['name'] = "Conformer A (reference)"
    mol2.info['name'] = "Conformer B (+15°)"
    mol3.info['name'] = "Conformer C (-15°)"

    # Rotate mol2/mol3 around origin (which is their COM after centering)
    def rot_z(atoms, angle):
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        atoms.set_positions(atoms.get_positions() @ R.T)

    rot_z(mol2, np.pi / 12)   # +15°
    rot_z(mol3, -np.pi / 12)  # -15°

    viewer = OverlayViewer([mol1, mol2, mol3], colorBy="Molecule")
    viewer.save_html(os.path.join(OUTPUT_DIR, "overlay_conformers.html"))
    print("Created overlay_conformers.html")


def create_overlay_colormap_viewer():
    """Create overlay with colormap."""
    # Create optimization-like trajectory (bond-length scan)
    water = molecule("H2O")
    # Center at origin so all frames share COM = [0,0,0]
    pos0 = water.get_positions()
    water.set_positions(pos0 - pos0.mean(axis=0))
    trajectory = []

    for i in range(10):
        atoms = water.copy()
        # Scale around origin (which is the COM after centering)
        atoms.set_positions(atoms.get_positions() * (1.0 + i * 0.02))
        trajectory.append(atoms)

    viewer = OverlayViewer(trajectory, colorBy="Colormap", colormap="viridis")
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
    """Create Au(111) + H2O adsorption relaxation trajectory."""
    # Create 3-layer Au(111) slab with 2x2 surface cell
    au_slab = fcc111('Au', size=(2, 2, 3), vacuum=8.0)

    # Add H2O molecule above the surface
    au_positions = au_slab.get_positions()
    z_max = au_positions[:, 2].max()

    # H2O initial position (high above surface)
    h2o_init_z = z_max + 4.0

    # Create relaxation trajectory (H2O approaching and relaxing on surface)
    trajectory = []
    n_frames = 20

    for i in range(n_frames):
        atoms = au_slab.copy()

        # H2O descends and relaxes
        progress = i / (n_frames - 1)
        # Exponential decay for smooth approach
        h2o_z = z_max + 2.2 + 1.8 * np.exp(-3 * progress)

        # Add H2O (O and two H atoms)
        center_xy = au_positions[:, :2].mean(axis=0)
        # Slight oscillation during relaxation
        wobble = 0.1 * np.sin(4 * np.pi * progress) * (1 - progress)

        h2o_positions = [
            [center_xy[0] + wobble, center_xy[1], h2o_z],  # O
            [center_xy[0] - 0.76 + wobble, center_xy[1] + 0.59, h2o_z + 0.2],  # H
            [center_xy[0] + 0.76 + wobble, center_xy[1] + 0.59, h2o_z + 0.2],  # H
        ]

        from ase import Atoms as AseAtoms
        h2o = AseAtoms('OH2', positions=h2o_positions)
        atoms += h2o

        # Energy: starts high, decreases as H2O adsorbs
        energy = -150.0 - 2.0 * progress + 0.3 * np.exp(-5 * progress) * np.sin(6 * np.pi * progress)
        atoms.info['energy'] = energy

        trajectory.append(atoms)

    viewer = MolecularViewer(trajectory, style="metallic", showCell=True, showEnergyPlot=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "gold_surface.html"))
    print("Created gold_surface.html (Au + H2O relaxation)")


def create_graphene_viewer():
    """Create graphene with phonon normal modes."""
    # Create small graphene nanoribbon
    graphene = graphene_nanoribbon(3, 3, type='zigzag', saturated=False, vacuum=5.0)
    n_atoms = len(graphene)

    # Create fake phonon modes (out-of-plane vibrations are characteristic)
    mode_vectors = []
    frequencies = []

    positions = graphene.get_positions()
    center = positions.mean(axis=0)

    # Mode 1: Breathing mode (in-plane radial)
    mode1 = []
    for pos in positions:
        r = pos - center
        r_norm = np.linalg.norm(r[:2])
        if r_norm > 0.1:
            disp = r[:2] / r_norm * 0.3
            mode1.append([disp[0], disp[1], 0.0])
        else:
            mode1.append([0.0, 0.0, 0.0])
    mode_vectors.append(mode1)
    frequencies.append(1580.0)  # G-band like

    # Mode 2: Out-of-plane acoustic (ZA mode)
    mode2 = []
    for pos in positions:
        # Wave-like displacement
        phase = 0.5 * (pos[0] + pos[1])
        mode2.append([0.0, 0.0, 0.4 * np.sin(phase)])
    mode_vectors.append(mode2)
    frequencies.append(450.0)

    # Mode 3: Optical out-of-plane (ZO mode)
    mode3 = []
    for i, pos in enumerate(positions):
        # Alternating up/down pattern
        sign = 1 if i % 2 == 0 else -1
        mode3.append([0.0, 0.0, 0.3 * sign])
    mode_vectors.append(mode3)
    frequencies.append(860.0)

    viewer = NormalViewer(
        graphene,
        mode_vectors=mode_vectors,
        frequencies=frequencies,
        showModeVector=True,
        style="cartoon"
    )
    viewer.save_html(os.path.join(OUTPUT_DIR, "graphene.html"))
    print("Created graphene.html (phonon modes)")


def create_nanotube_viewer():
    """Create carbon nanotube with vibrational modes."""
    cnt = nanotube(5, 0, length=2, vacuum=5.0)
    n_atoms = len(cnt)

    positions = cnt.get_positions()
    center = positions.mean(axis=0)

    mode_vectors = []
    frequencies = []

    # Mode 1: Radial breathing mode (RBM) - characteristic of CNTs
    mode1 = []
    for pos in positions:
        r = pos[:2] - center[:2]
        r_norm = np.linalg.norm(r)
        if r_norm > 0.1:
            disp = r / r_norm * 0.4
            mode1.append([disp[0], disp[1], 0.0])
        else:
            mode1.append([0.0, 0.0, 0.0])
    mode_vectors.append(mode1)
    frequencies.append(280.0)  # RBM frequency

    # Mode 2: Longitudinal acoustic mode
    mode2 = []
    z_min, z_max = positions[:, 2].min(), positions[:, 2].max()
    z_range = z_max - z_min if z_max > z_min else 1.0
    for pos in positions:
        z_rel = (pos[2] - z_min) / z_range
        mode2.append([0.0, 0.0, 0.3 * np.sin(np.pi * z_rel)])
    mode_vectors.append(mode2)
    frequencies.append(850.0)

    # Mode 3: G-band like tangential mode
    mode3 = []
    for pos in positions:
        r = pos[:2] - center[:2]
        r_norm = np.linalg.norm(r)
        if r_norm > 0.1:
            # Tangential direction
            tangent = np.array([-r[1], r[0]]) / r_norm * 0.25
            mode3.append([tangent[0], tangent[1], 0.0])
        else:
            mode3.append([0.0, 0.0, 0.0])
    mode_vectors.append(mode3)
    frequencies.append(1590.0)  # G-band

    viewer = NormalViewer(
        cnt,
        mode_vectors=mode_vectors,
        frequencies=frequencies,
        showModeVector=True,
        style="neon",
        backgroundColor="#000000"
    )
    viewer.save_html(os.path.join(OUTPUT_DIR, "carbon_nanotube.html"))
    print("Created carbon_nanotube.html (vibrational modes)")


def create_fcc_metal_viewer():
    """Create FCC copper crystal."""
    cu = bulk('Cu', 'fcc', a=3.61, cubic=True)
    cu = cu * (3, 3, 3)
    viewer = MolecularViewer(cu, style="metallic", showCell=True)
    viewer.save_html(os.path.join(OUTPUT_DIR, "copper_fcc.html"))
    print("Created copper_fcc.html")


def create_phosphine_charge_viewer():
    """Create triphenylphosphine molecule with partial charge visualization."""
    # Triphenylphosphine-like molecule with P center
    phosphine = Atoms(
        symbols=['C', 'C', 'O', 'P', 'H', 'C', 'C', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'H', 'H',
                 'C', 'C', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'H', 'H',
                 'C', 'C', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'H', 'H',
                 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'H'],
        positions=[
            [-1.450827, -1.737032, 1.036957],
            [-1.237239, -1.242182, 2.469770],
            [-0.302701, -0.154739, 2.173112],
            [-0.175123, -0.288858, 0.477351],
            [-0.692752, -1.973794, 3.079871],
            [-1.321641, 0.590023, -0.710356],
            [-1.416249, 0.316537, -2.082437],
            [-2.164868, 1.570102, -0.159671],
            [-2.333323, 1.007374, -2.883611],
            [-0.771664, -0.422915, -2.545414],
            [-3.086658, 2.252319, -0.957761],
            [-2.102498, 1.811603, 0.898579],
            [-3.175449, 1.972217, -2.325043],
            [-2.384614, 0.787877, -3.946900],
            [-3.731184, 3.004342, -0.510295],
            [-3.889341, 2.503663, -2.948197],
            [1.060058, 1.202925, 0.512177],
            [1.626613, 1.611358, -0.710060],
            [1.414893, 1.930588, 1.660761],
            [2.512089, 2.688546, -0.787579],
            [1.381869, 1.079259, -1.625809],
            [2.298839, 3.016923, 1.590660],
            [1.002143, 1.647801, 2.621096],
            [2.853216, 3.400310, 0.368173],
            [2.933792, 2.971811, -1.748789],
            [2.552737, 3.559290, 2.498335],
            [3.540268, 4.240706, 0.313683],
            [1.015325, -1.465638, -0.354058],
            [2.349228, -1.520924, 0.079158],
            [0.575887, -2.417359, -1.290946],
            [3.224214, -2.493662, -0.415966],
            [2.713772, -0.807609, 0.811543],
            [1.459158, -3.363106, -1.817026],
            [-0.462154, -2.437643, -1.607783],
            [2.786418, -3.408585, -1.376212],
            [4.248384, -2.527658, -0.054273],
            [1.104614, -4.075043, -2.557541],
            [3.468716, -4.155756, -1.772097],
            [-2.425796, -0.719459, 3.260493],
            [-2.086685, -0.280880, 4.205484],
            [-2.980400, 0.045097, 2.707630],
            [-3.113280, -1.539529, 3.498833],
            [-2.871566, -1.829498, 0.491716],
            [-3.449353, -2.593170, 1.031144],
            [-3.421248, -0.886713, 0.563562],
            [-2.862052, -2.121865, -0.564535],
            [-0.958633, -2.703662, 0.893846]
        ]
    )

    # Partial charges from DFT calculation
    charges = [
        0.0, 0.0, 0.06482878, 0.0, -0.01466073,
        0.0, 0.036884, 0.02864599, 0.04268293, -0.03649907,
        0.04313092, -0.01624414, 0.04600992, -0.01249285, -0.02825931,
        -0.01286479, 0.0, 0.03094519, 0.03708687, 0.04317401,
        -0.03103787, 0.04241059, -0.00997932, 0.04502739, -0.02730842,
        -0.02597865, -0.01794624, 0.01191878, 0.02973731, 0.03265756,
        0.0, -0.01961894, 0.04197893, -0.03086294, 0.04600057,
        -0.02808898, -0.02057052, -0.01400861, 0.0, -0.02920944,
        -0.0244752, -0.00137575, 0.01974016, -0.03187145, -0.01698247,
        -0.02150444, -0.00701329
    ]

    phosphine.arrays['charges'] = np.array(charges)

    viewer = MolecularViewer(phosphine, style="cartoon", colorBy="Charge")
    viewer.save_html(os.path.join(OUTPUT_DIR, "phosphine_charges.html"))
    print("Created phosphine_charges.html (partial charge visualization)")


def create_ethanol_charge_viewer():
    """Create ethanol molecule with partial charge visualization."""
    ethanol = molecule("CH3CH2OH")

    # Partial charges for ethanol (approximate)
    # Order: C(methyl), H, H, H, C(methylene), H, H, O, H
    charges = [
        -0.18,  # C (methyl)
        0.06,   # H
        0.06,   # H
        0.06,   # H
        0.15,   # C (connected to O)
        0.04,   # H
        0.04,   # H
        -0.65,  # O (hydroxyl)
        0.42,   # H (hydroxyl, acidic)
    ]

    ethanol.arrays['charges'] = np.array(charges)

    viewer = MolecularViewer(ethanol, style="glossy", colorBy="Charge")
    viewer.save_html(os.path.join(OUTPUT_DIR, "ethanol_charges.html"))
    print("Created ethanol_charges.html (partial charge visualization)")


def create_casio3_polyhedron_viewer():
    """Create CaSiO3 perovskite viewer with polyhedron visualization (inline structure)."""
    import io
    from ase.io import read

    # CaSiO3 tetragonal perovskite (I4/mcm) - generated using pymatgen
    cif_content = """\
# generated using pymatgen
data_CaSiO3
_symmetry_space_group_name_H-M   I4/mcm
_cell_length_a   5.07858400
_cell_length_b   5.07858400
_cell_length_c   7.28534200
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   140
_chemical_formula_structural   CaSiO3
_chemical_formula_sum   'Ca4 Si4 O12'
_cell_volume   187.90365339
_cell_formula_units_Z   4
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
  2  '-x, -y, -z'
  3  '-y, x, z'
  4  'y, -x, -z'
  5  '-x, -y, z'
  6  'x, y, -z'
  7  'y, -x, z'
  8  '-y, x, -z'
  9  'x, -y, -z+1/2'
  10  '-x, y, z+1/2'
  11  '-y, -x, -z+1/2'
  12  'y, x, z+1/2'
  13  '-x, y, -z+1/2'
  14  'x, -y, z+1/2'
  15  'y, x, -z+1/2'
  16  '-y, -x, z+1/2'
  17  'x+1/2, y+1/2, z+1/2'
  18  '-x+1/2, -y+1/2, -z+1/2'
  19  '-y+1/2, x+1/2, z+1/2'
  20  'y+1/2, -x+1/2, -z+1/2'
  21  '-x+1/2, -y+1/2, z+1/2'
  22  'x+1/2, y+1/2, -z+1/2'
  23  'y+1/2, -x+1/2, z+1/2'
  24  '-y+1/2, x+1/2, -z+1/2'
  25  'x+1/2, -y+1/2, -z'
  26  '-x+1/2, y+1/2, z'
  27  '-y+1/2, -x+1/2, -z'
  28  'y+1/2, x+1/2, z'
  29  '-x+1/2, y+1/2, -z'
  30  'x+1/2, -y+1/2, z'
  31  'y+1/2, x+1/2, -z'
  32  '-y+1/2, -x+1/2, z'
loop_
 _atom_type_symbol
 _atom_type_oxidation_number
  Ca2+  2.0
  Si4+  4.0
  O2-  -2.0
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Ca2+  Ca0  4  0.00000000  0.50000000  0.25000000  1
  Si4+  Si1  4  0.00000000  0.00000000  0.00000000  1
  O2-  O2  8  0.22204050  0.27795950  0.00000000  1
  O2-  O3  4  0.00000000  0.00000000  0.25000000  1
"""
    import tempfile, os
    with tempfile.NamedTemporaryFile(suffix='.cif', mode='w', delete=False) as f:
        f.write(cif_content)
        tmp_path = f.name
    try:
        atoms = read(tmp_path)
    finally:
        os.unlink(tmp_path)

    atoms = atoms * (2, 2, 2)
    viewer = MolecularViewer(
        atoms,
        style="default",
        showCell=True,
        showBond=True,
        showPolyhedron=True,
        polyhedronOpacity=0.25,
        atomSize=0.35,
        bondThreshold=1.1
    )
    viewer.save_html(os.path.join(OUTPUT_DIR, "casio3_polyhedron.html"))
    print("Created casio3_polyhedron.html")


def create_frag_selector_basic():
    """Create basic FragSelector example (methanol)."""
    atoms = molecule('CH3OH')
    viewer = FragSelector(atoms, style='cartoon')
    viewer.save_html(os.path.join(OUTPUT_DIR, "frag_selector_basic.html"))
    print("Created frag_selector_basic.html")


def create_frag_selector_multifrag():
    """Create multi-fragment FragSelector example (two water molecules)."""
    water1 = molecule('H2O')
    water2 = molecule('H2O')
    water2.translate([5, 0, 0])
    system = water1 + water2
    viewer = FragSelector(system, style='cartoon')
    viewer.save_html(os.path.join(OUTPUT_DIR, "frag_selector_multifrag.html"))
    print("Created frag_selector_multifrag.html")


def create_frag_selector_bis_iodo():
    """Create large conjugated system FragSelector example (bis-iodo benzimidazole dimer)."""
    from ase.io import read
    xyz_path = os.path.join(os.path.dirname(__file__), "assets", "examples", "dimer.xyz")
    atoms = read(xyz_path)
    viewer = FragSelector(atoms, style='cartoon')
    viewer.save_html(os.path.join(OUTPUT_DIR, "frag_selector_bis_iodo.html"))
    print("Created frag_selector_bis_iodo.html")


def create_carbazole_ring_viewer():
    """Create carbazole molecule viewer with ring highlighting (inline coordinates)."""
    import os
    from ase import Atoms

    # Carbazole (C12H9N core + phosphonate group) - 45 atoms
    symbols = [
        'C','C','C','C','C','C','C','N','C','C','C','C','C','C','C',
        'C','C','C','C','P','O','O','O',
        'H','H','H','H','H','H','H','H','H','H','H','H',
        'H','H','H','H','H','H','H','H','H','H',
    ]
    positions = [
        [-1.268587,  3.133147, 10.747185],
        [-0.910866,  2.368808,  9.506540],
        [-0.690901,  0.983936,  9.572118],
        [-0.384794,  0.255745,  8.401030],
        [-0.315087,  0.936147,  7.185078],
        [-0.550664,  2.317767,  7.095945],
        [-0.852236,  3.015722,  8.270608],
        [ 0.000000,  0.000000,  6.218598],
        [ 0.118345, -1.260468,  6.773151],
        [-0.115645, -1.124745,  8.141747],
        [-0.052302, -2.274078,  8.960522],
        [ 0.226774, -3.529192,  8.397141],
        [ 0.466183, -3.628371,  7.025164],
        [ 0.410736, -2.505323,  6.192375],
        [ 0.325191, -4.742434,  9.274520],
        [ 0.188517,  0.300979,  4.805273],
        [-1.114221,  0.173663,  4.016158],
        [-0.970087,  0.574017,  2.547094],
        [-0.091333, -0.375951,  1.738248],
        [ 0.000000,  0.000000,  0.000000],
        [ 0.846897, -0.876297, -0.857630],
        [-1.510074,  0.075719, -0.524504],
        [ 0.439665,  1.536509, -0.097183],
        [-0.797589,  2.687782, 11.629830],
        [-2.353966,  3.132756, 10.886772],
        [-0.920085,  4.169274, 10.682103],
        [-0.754857,  0.462445, 10.523085],
        [-0.508846,  2.843024,  6.147438],
        [-1.042825,  4.085461,  8.207878],
        [-0.224754, -2.178096, 10.028901],
        [ 0.700456, -4.595056,  6.583224],
        [ 0.592884, -2.609213,  5.127879],
        [-0.346025, -4.657379, 10.135616],
        [ 1.351415, -4.860402,  9.635738],
        [ 0.036392, -5.645590,  8.726657],
        [ 0.598193,  1.313979,  4.727826],
        [ 0.958922, -0.377103,  4.424704],
        [-1.878057,  0.815466,  4.473496],
        [-1.497217, -0.852179,  4.083605],
        [-0.574534,  1.593850,  2.475698],
        [-1.968957,  0.600580,  2.095905],
        [ 0.930434, -0.364498,  2.131032],
        [-0.465076, -1.400279,  1.846490],
        [-1.576519,  0.092936, -1.497010],
        [ 0.698440,  1.801183, -0.999255],
    ]
    atoms = Atoms(symbols=symbols, positions=positions)
    viewer = MolecularViewer(
        atoms,
        style="cartoon",
        showRings=True,
        ringOpacity=0.35,
        showBond=True,
        bondThreshold=1.15
    )
    viewer.save_html(os.path.join(OUTPUT_DIR, "carbazole_rings.html"))
    print("Created carbazole_rings.html")


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
    # Charge visualization
    create_phosphine_charge_viewer()
    create_ethanol_charge_viewer()
    # Polyhedron & Ring highlight
    create_casio3_polyhedron_viewer()
    create_carbazole_ring_viewer()
    # Fragment selector
    create_frag_selector_basic()
    create_frag_selector_multifrag()
    create_frag_selector_bis_iodo()
    print("Done!")
