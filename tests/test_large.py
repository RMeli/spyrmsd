from spyrmsd import io, spyrmsd

import requests
import os
import zipfile
import warnings

import numpy as np

import pytest

n_systems = 4554

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
    ofname = os.path.join(tlpath, "docking.zip") # Downloaded .zip file
    odname = os.path.join(tlpath, "docking") # Data directory

    def unzip():
        with zipfile.ZipFile(ofname,"r") as zipf:
            zipf.extractall(os.path.join(tlpath, "docking"))

    if os.path.isdir(odname):
        pass # Folder already exists, nothing to do
    elif os.path.isfile(ofname):
        unzip()
    else:
        # Archive on Zenodo
        version = "0.1.0"
        url="https://github.com/RMeli/PDBbind-docking/archive/master.zip"

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
    assert len(download) == n_systems

    for id in download:
        assert os.path.isdir(path_from_id(id, path))

@pytest.mark.large
@pytest.mark.parametrize("minimize", [True, False])
@pytest.mark.parametrize(
    "idx", np.random.randint(0, n_systems, size=250)
)
def test_rmsd(idx, download, path, minimize):
    id = download[idx]

    p = path_from_id(id, path)

    results = np.loadtxt(os.path.join(p, "obrms-min.dat" if minimize else "obrms.dat"))

    ref = io.loadmol(os.path.join(p, f"{id}_ligand.sdf"))

    try:
        mols = io.loadallmols(os.path.join(p, f"{id}_dock.sdf"))
    except OSError: # Docking didn't find any configuration for some systems
        warnings.warn(f"File {id}_dock.sdf not found.", RuntimeWarning)
        return

    rmsds = spyrmsd.rmsdwrapper(ref, mols, minimize=minimize, strip=True)

    assert np.allclose(rmsds, results)