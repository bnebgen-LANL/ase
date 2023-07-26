"""Test write and read."""
import io
import re

import pytest
from ase import Atoms
from ase.build import bulk
from ase.data import atomic_numbers
from ase.io.lammpsdata import read_lammps_data, write_lammps_data


def compare(atoms: Atoms, atoms_ref: Atoms):
    """Compare two `Atoms` objects"""
    assert all(atoms.numbers == atoms_ref.numbers)
    assert atoms.get_masses() == pytest.approx(atoms_ref.get_masses())

    # Note: Raw positions cannot be compared.
    # `write_lammps_data` changes cell orientation.
    assert atoms.get_scaled_positions() == pytest.approx(
        atoms_ref.get_scaled_positions())


@pytest.mark.parametrize("masses", [False, True])
class _Base:
    def run(self, atoms_ref: Atoms, masses: bool):
        """Run tests"""
        self.check_explicit_numbers(atoms_ref, masses)
        if masses:
            self.check_masses2numbers(atoms_ref)

    def check_explicit_numbers(self, atoms_ref: Atoms, masses: bool):
        """Check if write-read is consistent when giving Z_of_type."""
        buf = io.StringIO()
        write_lammps_data(buf, atoms_ref, masses=masses)
        buf.seek(0)
        # By default, write_lammps_data assigns atom types to the elements in
        # alphabetical order. To be consistent, here spiecies are also sorted.
        species = sorted(set(atoms_ref.get_chemical_symbols()))
        Z_of_type = {i + 1: atomic_numbers[s] for i, s in enumerate(species)}
        atoms = read_lammps_data(buf, Z_of_type=Z_of_type, style="atomic")
        compare(atoms, atoms_ref)

    def check_masses2numbers(self, atoms_ref: Atoms):
        """Check if write-read is consistent when guessing atomic numbers."""
        buf = io.StringIO()
        write_lammps_data(buf, atoms_ref, masses=True)
        buf.seek(0)
        atoms = read_lammps_data(buf, style="atomic")
        compare(atoms, atoms_ref)


@pytest.mark.parametrize("cubic", [False, True])
class TestCubic(_Base):
    """Test cubic structures."""

    def test_bcc(self, cubic: bool, masses: bool):
        """Test bcc."""
        atoms_ref = bulk("Li", "bcc", cubic=cubic)
        self.run(atoms_ref, masses)

    def test_fcc(self, cubic: bool, masses: bool):
        """Test fcc."""
        atoms_ref = bulk("Cu", "fcc", cubic=cubic)
        self.run(atoms_ref, masses)

    def test_rocksalt(self, cubic: bool, masses: bool):
        """Test rocksalt."""
        atoms_ref = bulk("NaCl", "rocksalt", a=1.0, cubic=cubic)
        self.run(atoms_ref, masses)

    def test_fluorite(self, cubic: bool, masses: bool):
        """Test fluorite."""
        atoms_ref = bulk("CaF2", "fluorite", a=1.0, cubic=cubic)
        self.run(atoms_ref, masses)


@pytest.mark.parametrize("orthorhombic", [False, True])
class TestOrthorhombic(_Base):
    """Test orthorhombic structures."""

    def test_hcp(self, orthorhombic: bool, masses: bool):
        """Test hcp."""
        atoms_ref = bulk("Mg", "hcp", orthorhombic=orthorhombic)
        self.run(atoms_ref, masses)


@pytest.mark.parametrize("atom_style", ["atomic", "charge", "full"])
def test_atom_style(atom_style: str):
    """Test `atom_style`"""
    atoms_ref = bulk("Cu", "fcc")

    # note that we should do `masses=True` to later guess atomic numbers
    buf = io.StringIO()
    write_lammps_data(buf, atoms_ref, masses=True, atom_style=atom_style)

    # test if the `atom_style` can be guessed from the comment
    buf.seek(0)
    atoms = read_lammps_data(buf, style=None)
    compare(atoms, atoms_ref)

    # test when `atom_style` is explicitly specified
    buf.seek(0)
    atoms = read_lammps_data(buf, style=atom_style)
    compare(atoms, atoms_ref)

    if atom_style not in ["atomic", "full"]:
        return

    # test if `atom_style` can be guessed from the length of fields
    # remove comment in the "Atoms" line
    buf = io.StringIO(re.sub(".*Atoms.*#.*", "Atoms", buf.getvalue()))
    atoms = read_lammps_data(buf, style=None)
    compare(atoms, atoms_ref)
