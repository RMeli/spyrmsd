from spyrmsd import io, rmsd

import pytest

import numpy as np
import os

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, "data", "docking")


def systems():
    return os.listdir(molpath)


@pytest.fixture(scope="module", params=systems())
def molecules(request):
    system = request.param

    fref = os.path.join(molpath, system, f"{system}_ligand.sdf")
    fmols = os.path.join(molpath, system, f"{system}_dock.sdf")

    ref = io.load(os.path.join(molpath, fref))
    mols = io.loadall(os.path.join(molpath, fmols))

    ref = io.to_molecule(ref)
    mols = [io.to_molecule(mol) for mol in mols]

    ref.strip()
    for mol in mols:
        mol.strip()

    return ref, mols, system


@pytest.mark.parametrize("cache", [True, False])
def test_benchmark_symmrmsd(cache, molecules, benchmark):

    ref, mols, system = molecules

    coords = [mol.coordinates for mol in mols]

    RMSDs = benchmark(
        rmsd.multirmsd_isomorphic,
        ref.coordinates,
        coords,
        ref.adjacency_matrix,
        mols[0].adjacency_matrix,
        ref.atomicnums,
        mols[0].atomicnums,
        cache=cache,
    )

    testRMSDs = np.loadtxt(os.path.join(molpath, system, "obrms.dat"))

    assert np.allclose(RMSDs, testRMSDs)
