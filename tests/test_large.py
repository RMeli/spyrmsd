import os
import zipfile

import numpy as np
import pytest
import requests

from spyrmsd import io, qcp, rmsd


@pytest.fixture(autouse=True, params=[True, False])
def lambda_max_failure(monkeypatch, request):
    """
    Monkey patch fixture for :code:`lambda_max` function to simulate convergence
    failures.

    Notes
    -----
    The :fun:`lambda_max` function can sometimes fail to converge and raises a
    :code:`RuntimeError` (see GitHub Issue #35 by @kjelljorner). If this occours,
    there is an automatic fallback to the explicit calculation of the highest
    eigenvalue. This is not easy to reproduce but need to be tested.

    Using monkey patching we run all tests with the original :fun:`lambda_max` function
    and again with a patched :fun:`lambda_max` function that always raise the
    :code:`RuntimeError`, so that the fallback is tested instead.

    https://github.com/RMeli/spyrmsd/issues/35
    """
    if request.param:
        # Patch lambda_max to always raise an exception
        # This enforces _lambda_max_eig to be used instead
        def lambda_max_failure(Ga, Gb, c2, c1, c0):
            # Simulate Newton method convergence failure
            raise RuntimeError

        monkeypatch.setattr(qcp, "lambda_max", lambda_max_failure)


@pytest.fixture
def tlpath():
    """
    test_long.py path.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def path(tlpath):
    """
    Path to docking files.
    """
    return os.path.join(tlpath, "docking", "PDBbind-docking-master", "docking")


@pytest.fixture
def download(tlpath, path):
    ofname = os.path.join(tlpath, "docking.zip")  # Downloaded .zip file
    odname = os.path.join(tlpath, "docking")  # Data directory

    def unzip():
        with zipfile.ZipFile(ofname, "r") as zipf:
            zipf.extractall(os.path.join(tlpath, "docking"))

    if os.path.isdir(odname):
        pass  # Folder already exists, nothing to do
    elif os.path.isfile(ofname):
        unzip()
    else:
        # TODO: Change to Zenodo archive
        url = "https://github.com/RMeli/PDBbind-docking/archive/master.zip"

        r = requests.get(url)

        with open(ofname, "wb") as fout:
            fout.write(r.content)

        unzip()

    # Return list of systems
    return [id for id in os.listdir(path) if len(id) == 4]


def path_from_id(id: str, path: str):
    return os.path.join(path, id)


@pytest.mark.large
def test_dowload(download, path):
    assert os.path.isdir(path)
    assert len(download) == pytest.n_systems

    for id in download:
        assert os.path.isdir(path_from_id(id, path))


@pytest.mark.large
@pytest.mark.parametrize("minimize", [True, False])
def test_rmsd(idx, download, path, minimize):
    id = download[idx]

    p = path_from_id(id, path)

    try:
        ref = io.loadmol(os.path.join(p, f"{id}_ligand.sdf"))
        mols = io.loadallmols(os.path.join(p, f"{id}_dock.sdf"))

        # Load results obtained with OpenBabel
        results = np.loadtxt(
            os.path.join(p, "obrms-min.dat" if minimize else "obrms.dat")
        )

    except OSError:  # Docking didn't find any configuration for some systems
        pytest.xfail("Problems loading PDB ID {id}.")

    # Some molecules in the dataset give errors when parsed with RDKit
    if ref is None or None in mols:
        pytest.xfail("Problems loading PDB ID {id}.")

    rmsds = rmsd.rmsdwrapper(ref, mols, minimize=minimize, strip=True)

    assert np.allclose(rmsds, results)
