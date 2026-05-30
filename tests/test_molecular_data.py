"""Test MolecularData conversion between ASE Atoms and dict."""

import pytest
import numpy as np
from ase import Atoms

from aseview.wrapper import MolecularData


@pytest.fixture
def water():
    """Simple water molecule."""
    return Atoms(
        "H2O",
        positions=[[0.0, 0.0, 0.0], [0.0, 0.76, 0.59], [0.0, -0.76, 0.59]],
    )


@pytest.fixture
def bulk_si():
    """Silicon crystal with periodic boundary conditions."""
    from ase.build import bulk
    return bulk("Si", "diamond", a=5.43)


class TestMolecularDataFromAtoms:
    def test_basic_conversion(self, water):
        data = MolecularData.from_atoms(water)
        assert "positions" in data
        assert "symbols" in data
        assert len(data["positions"]) == 3
        assert data["symbols"] == ["H", "H", "O"]

    def test_positions_are_lists(self, water):
        data = MolecularData.from_atoms(water)
        assert isinstance(data["positions"], list)
        assert isinstance(data["positions"][0], list)

    def test_cell_included_for_periodic(self, bulk_si):
        data = MolecularData.from_atoms(bulk_si)
        assert "cell" in data
        assert len(data["cell"]) == 3

    def test_no_cell_for_molecule(self, water):
        data = MolecularData.from_atoms(water)
        assert "cell" not in data

    def test_charges_preserved(self, water):
        charges = np.array([0.4, -0.2, -0.2])
        water.arrays["charges"] = charges
        data = MolecularData.from_atoms(water)
        assert "charges" in data
        assert len(data["charges"]) == 3
        assert abs(data["charges"][0] - 0.4) < 1e-6

    def test_magmoms_preserved_from_info(self, water):
        water.info["magmom"] = [1.0, -1.0, 0.0]
        data = MolecularData.from_atoms(water)
        assert "magmoms" in data
        assert len(data["magmoms"]) == 3
        assert abs(data["magmoms"][1] + 1.0) < 1e-6

    def test_magmoms_preserved_from_initial_array(self, water):
        water.set_initial_magnetic_moments([0.5, -0.5, 0.0])
        data = MolecularData.from_atoms(water)
        assert "magmoms" in data
        assert len(data["magmoms"]) == 3
        assert abs(data["magmoms"][0] - 0.5) < 1e-6

    def test_vector_magmoms_are_converted_to_magnitudes(self, water):
        water.info["magmoms"] = [
            [1.0, 0.0, 0.0],
            [0.0, -2.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
        data = MolecularData.from_atoms(water)
        assert "magmoms" in data
        np.testing.assert_allclose(data["magmoms"], [1.0, 2.0, 0.0])


class TestMolecularDataToAtoms:
    def test_roundtrip(self, water):
        data = MolecularData.from_atoms(water)
        atoms = MolecularData.to_atoms(data)
        assert len(atoms) == 3
        assert atoms.get_chemical_symbols() == ["H", "H", "O"]
        np.testing.assert_allclose(
            atoms.get_positions(),
            water.get_positions(),
            atol=1e-10,
        )

    def test_roundtrip_periodic(self, bulk_si):
        data = MolecularData.from_atoms(bulk_si)
        atoms = MolecularData.to_atoms(data)
        assert atoms.pbc.all()
        np.testing.assert_allclose(
            atoms.get_cell()[:],
            bulk_si.get_cell()[:],
            atol=1e-10,
        )
