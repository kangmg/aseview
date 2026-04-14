"""
Hessian file parsers for various quantum chemistry programs.

Supported formats:
- ORCA: .hess files
- VASP: OUTCAR files (IBRION=5 or 6)
"""

import re
import numpy as np
from typing import Tuple, List, Optional
from pathlib import Path


def parse_orca_hess(filepath: str) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Parse ORCA .hess file to extract normal modes and frequencies.

    Args:
        filepath: Path to the ORCA .hess file

    Returns:
        Tuple of (frequencies, normal_modes, n_atoms)
        - frequencies: 1D array of vibrational frequencies in cm^-1
        - normal_modes: 2D array of shape (n_modes, 3*n_atoms), each row is a mode vector
        - n_atoms: Number of atoms
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Hessian file not found: {filepath}")

    with open(filepath, "r") as f:
        content = f.read()

    # Parse $vibrational_frequencies section
    frequencies = _parse_orca_frequencies(content)

    # Parse $normal_modes section
    normal_modes = _parse_orca_normal_modes(content)

    # Calculate number of atoms
    n_coords = normal_modes.shape[1]
    n_atoms = n_coords // 3

    return frequencies, normal_modes, n_atoms


def _parse_orca_frequencies(content: str) -> np.ndarray:
    """Parse the $vibrational_frequencies section."""
    # Find the section
    start_marker = "$vibrational_frequencies"
    end_markers = ["$", "#"]

    start_idx = content.find(start_marker)
    if start_idx == -1:
        raise ValueError("Could not find $vibrational_frequencies section in ORCA .hess file")

    # Find the section content
    section_start = start_idx + len(start_marker)
    section_end = len(content)

    for marker in end_markers:
        next_section = content.find(marker, section_start + 1)
        if next_section != -1 and next_section < section_end:
            section_end = next_section

    section = content[section_start:section_end].strip()
    lines = section.split("\n")

    # First line is the number of frequencies
    n_freqs = int(lines[0].strip())

    frequencies = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            # Format: index  frequency
            freq = float(parts[1])
            frequencies.append(freq)

    if len(frequencies) != n_freqs:
        raise ValueError(f"Expected {n_freqs} frequencies, got {len(frequencies)}")

    return np.array(frequencies)


def _parse_orca_normal_modes(content: str) -> np.ndarray:
    """Parse the $normal_modes section."""
    start_marker = "$normal_modes"
    end_markers = ["$", "#"]

    start_idx = content.find(start_marker)
    if start_idx == -1:
        raise ValueError("Could not find $normal_modes section in ORCA .hess file")

    # Find the section content
    section_start = start_idx + len(start_marker)
    section_end = len(content)

    for marker in end_markers:
        next_section = content.find(marker, section_start + 1)
        if next_section != -1 and next_section < section_end:
            section_end = next_section

    section = content[section_start:section_end].strip()
    lines = section.split("\n")

    # First line contains dimensions: n_coords n_modes
    dims = lines[0].strip().split()
    n_coords = int(dims[0])
    n_modes = int(dims[1])

    # Initialize the normal modes matrix
    normal_modes = np.zeros((n_modes, n_coords))

    # Parse the mode data (block format)
    # ORCA outputs modes in blocks of up to 5 columns
    current_line = 1
    mode_col_start = 0

    while mode_col_start < n_modes:
        # Skip empty lines
        while current_line < len(lines) and not lines[current_line].strip():
            current_line += 1

        if current_line >= len(lines):
            break

        # Header line with mode indices
        header = lines[current_line].strip().split()
        n_cols_in_block = len(header)
        current_line += 1

        # Read coordinate rows
        for coord_idx in range(n_coords):
            if current_line >= len(lines):
                break

            line = lines[current_line].strip()
            if not line:
                current_line += 1
                continue

            parts = line.split()
            # First part is coordinate index, rest are mode values
            for col_offset, value in enumerate(parts[1 : 1 + n_cols_in_block]):
                mode_idx = mode_col_start + col_offset
                if mode_idx < n_modes:
                    normal_modes[mode_idx, coord_idx] = float(value)

            current_line += 1

        mode_col_start += n_cols_in_block

    return normal_modes


def parse_vasp_outcar(filepath: str) -> Tuple[np.ndarray, np.ndarray, int, np.ndarray]:
    """
    Parse VASP OUTCAR file to extract normal modes and frequencies.

    The OUTCAR must contain the "Eigenvectors and eigenvalues of the dynamical
    matrix" section produced by IBRION=5 or IBRION=6 calculations.

    For calculations with selective dynamics (only a subset of atoms free),
    the returned n_atoms and free_atom_positions reflect only the free atoms.
    Use the free_atom_positions to match these atoms to a full structure.

    Args:
        filepath: Path to the VASP OUTCAR file

    Returns:
        Tuple of (frequencies, normal_modes, n_atoms, free_atom_positions)
        - frequencies: 1D array of vibrational frequencies in cm^-1
          (negative values indicate imaginary / soft modes)
        - normal_modes: 2D array of shape (n_modes, 3*n_atoms)
        - n_atoms: Number of free atoms in the eigenvector block
        - free_atom_positions: 2D array of shape (n_atoms, 3) — Cartesian
          coordinates of the free atoms as listed in the eigenvector block
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"OUTCAR file not found: {filepath}")

    with open(filepath, "r") as f:
        content = f.read()

    marker = "Eigenvectors and eigenvalues of the dynamical matrix"
    # Use rfind so we always get the final (converged) section
    idx = content.rfind(marker)
    if idx == -1:
        raise ValueError(
            'Could not find "Eigenvectors and eigenvalues of the dynamical matrix" '
            "in OUTCAR. Make sure IBRION=5 or IBRION=6 was used."
        )

    lines = content[idx:].split("\n")

    # Regex for mode-header lines:
    #   "   1 f  =   93.265082 THz ..."
    #   "   1 f/i=    2.714783 THz ..."
    mode_header_re = re.compile(
        r"^\s*\d+\s+f(/i)?\s*=\s*[\d.]+\s+THz\s+[\d.]+\s+2PiTHz\s+([\d.]+)\s+cm-1"
    )
    # Atom-displacement lines: 6 floats  (X  Y  Z  dx  dy  dz)
    atom_line_re = re.compile(
        r"^\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$"
    )

    frequencies: List[float] = []
    all_modes: List[List[List[float]]] = []
    all_positions: List[List[List[float]]] = []   # XYZ per mode
    current_disps: List[List[float]] = []
    current_pos: List[List[float]] = []
    in_mode = False

    for line in lines:
        m = mode_header_re.match(line)
        if m:
            # Save the previous mode
            if in_mode and current_disps:
                all_modes.append(current_disps)
                all_positions.append(current_pos)
                current_disps = []
                current_pos = []
            imaginary = m.group(1) is not None   # "/i" present
            freq = float(m.group(2))
            frequencies.append(-freq if imaginary else freq)
            in_mode = True
            continue

        if in_mode:
            am = atom_line_re.match(line)
            if am:
                x, y, z = float(am.group(1)), float(am.group(2)), float(am.group(3))
                dx, dy, dz = float(am.group(4)), float(am.group(5)), float(am.group(6))
                current_pos.append([x, y, z])
                current_disps.append([dx, dy, dz])

    # Flush the last mode
    if in_mode and current_disps:
        all_modes.append(current_disps)
        all_positions.append(current_pos)

    if not all_modes:
        raise ValueError(
            "No normal modes found in OUTCAR. "
            "Make sure IBRION=5 or IBRION=6 was used and the calculation completed."
        )

    n_atoms = len(all_modes[0])
    n_modes = len(all_modes)

    normal_modes = np.zeros((n_modes, n_atoms * 3))
    for mode_idx, mode in enumerate(all_modes):
        for atom_idx, (dx, dy, dz) in enumerate(mode[:n_atoms]):
            normal_modes[mode_idx, atom_idx * 3]     = dx
            normal_modes[mode_idx, atom_idx * 3 + 1] = dy
            normal_modes[mode_idx, atom_idx * 3 + 2] = dz

    # Use positions from the first mode as the canonical free-atom positions
    free_positions = np.array(all_positions[0]) if all_positions else np.zeros((n_atoms, 3))

    return np.array(frequencies), normal_modes, n_atoms, free_positions


def detect_hessian_format(filepath: str) -> str:
    """
    Auto-detect hessian file format from filename and content.

    Returns:
        'orca' or 'vasp'
    """
    name = Path(filepath).name.upper()
    if "OUTCAR" in name:
        return "vasp"

    # Fall back to content sniffing (read first 8 KB)
    with open(filepath, "r", errors="replace") as f:
        head = f.read(8192)

    if "Eigenvectors and eigenvalues of the dynamical matrix" in head:
        return "vasp"
    if "$vibrational_frequencies" in head:
        return "orca"

    raise ValueError(
        f"Cannot determine hessian format for '{filepath}'. "
        "Expected ORCA .hess or VASP OUTCAR."
    )


def reshape_modes_to_atoms(normal_modes: np.ndarray, n_atoms: int) -> List[List[List[float]]]:
    """
    Reshape flat normal modes to per-atom displacement vectors.

    Args:
        normal_modes: 2D array of shape (n_modes, 3*n_atoms)
        n_atoms: Number of atoms

    Returns:
        List of modes, where each mode is a list of [x, y, z] displacements per atom
        Shape: (n_modes, n_atoms, 3)
    """
    n_modes = normal_modes.shape[0]
    reshaped = []

    for mode_idx in range(n_modes):
        mode = normal_modes[mode_idx]
        atom_displacements = []
        for atom_idx in range(n_atoms):
            x = mode[3 * atom_idx]
            y = mode[3 * atom_idx + 1]
            z = mode[3 * atom_idx + 2]
            atom_displacements.append([x, y, z])
        reshaped.append(atom_displacements)

    return reshaped


def get_real_vibrations(
    frequencies: np.ndarray, normal_modes: np.ndarray, skip_translations_rotations: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Filter out translational and rotational modes (first 5-6 modes with ~0 frequency).

    Args:
        frequencies: Array of frequencies
        normal_modes: Array of normal modes
        skip_translations_rotations: If True, skip modes with frequency < 10 cm^-1

    Returns:
        Filtered (frequencies, normal_modes)
    """
    if not skip_translations_rotations:
        return frequencies, normal_modes

    # Find indices of real vibrations (frequency > threshold)
    threshold = 10.0  # cm^-1
    real_indices = np.where(np.abs(frequencies) > threshold)[0]

    return frequencies[real_indices], normal_modes[real_indices]
