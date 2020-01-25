import os

import pytest

from spyrmsd import io

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, "data/molecules/")


@pytest.mark.parametrize(
    "molfile, natoms, nbonds",
    [("benzene.sdf", 12, 12), ("ethanol.sdf", 9, 8), ("dialanine.sdf", 23, 22)],
)
def test_load_sdf(molfile, natoms: int, nbonds: int) -> None:

    m = io.load(os.path.join(molpath, molfile))

    assert io.numatoms(m) == natoms
    assert io.numbonds(m) == nbonds


@pytest.mark.parametrize(
    "molfile, natoms, nbonds", [("1cbr_ligand.mol2", 49, 49)],
)
def test_load_mol2(molfile, natoms: int, nbonds: int) -> None:

    m = io.load(os.path.join(molpath, molfile))

    assert io.numatoms(m) == natoms
    assert io.numbonds(m) == nbonds


@pytest.mark.parametrize(
    "molfile, natoms, nbonds",
    [("trp0.pdb", 217, 224), ("trp1.pdb", 217, 224), ("trp2.pdb", 217, 224)],
)
def test_load_pdb(molfile, natoms: int, nbonds: int) -> None:

    m = io.load(os.path.join(molpath, molfile))

    assert io.numatoms(m) == natoms
    assert io.numbonds(m) == nbonds


@pytest.mark.parametrize(
    "molfile, natoms, nbonds", [("1cbr_docking.sdf", 22, 22)],
)
def test_loadall_sdf(molfile, natoms: int, nbonds: int) -> None:

    ms = io.loadall(os.path.join(molpath, molfile))

    assert len(ms) == 10

    for m in ms:
        assert io.numatoms(m) == natoms
        assert io.numbonds(m) == nbonds


@pytest.mark.parametrize(
    "molfile, natoms, nbonds", [("1cbr_docking.mol2", 22, 22)],
)
def test_loadall_mol2(molfile, natoms: int, nbonds: int) -> None:

    try:
        ms = io.loadall(os.path.join(molpath, molfile))

        assert len(ms) == 10

        for m in ms:
            assert io.numatoms(m) == natoms
            assert io.numbonds(m) == nbonds
    except NotImplementedError:  # Mol2MolSupplier in RDkit is not supported
        pass  # TODO: Warning


@pytest.mark.parametrize(
    "molfile, natoms, nbonds", [("1cbr_docking.pdb", 22, 22)],
)
def test_loadall_pdb(molfile, natoms: int, nbonds: int) -> None:

    try:
        ms = io.loadall(os.path.join(molpath, molfile))

        assert len(ms) == 10

        for m in ms:
            assert io.numatoms(m) == natoms
            assert io.numbonds(m) == nbonds

    except NotImplementedError:  # PDBMolSupplier in RDkit is not supported
        pass  # Warning
