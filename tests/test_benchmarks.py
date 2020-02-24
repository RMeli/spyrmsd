import os
import re

import numpy as np
import pytest

from spyrmsd import io, rmsd

fdir = os.path.dirname(os.path.abspath(__file__))
molpath = os.path.join(fdir, "data", "docking")


def systems():
    dirs = os.listdir(molpath)

    # Match 4-character strings starting with a number
    systems = [d for d in dirs if re.match("[1-9]...", d)]

    return systems


@pytest.mark.slow
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


@pytest.mark.slow
@pytest.mark.parametrize("cache", [True, False])
def test_benchmark_symmrmsd(cache, molecules, benchmark):

    ref, mols, system = molecules

    coords = [mol.coordinates for mol in mols]

    RMSDs = benchmark(
        rmsd.symmrmsd,
        ref.coordinates,
        coords,
        ref.atomicnums,
        mols[0].atomicnums,
        ref.adjacency_matrix,
        mols[0].adjacency_matrix,
        cache=cache,
    )

    testRMSDs = np.loadtxt(os.path.join(molpath, system, "obrms.dat"))

    assert np.allclose(RMSDs, testRMSDs)
