from pyrmsd import utils

import pytest
import numpy as np

def test_format():
    for ext in ["pdf", "smi", "xyz"]:
        for fname in ["test", "root/test", "root/test.test"]:
            fmt = utils.format(f"{fname}.{ext}")
            assert fmt == ext

def test_format_openbabel():

    for ext in ["smi", "pdb"]:
        for fname in ["test", "root/test", "root/test.test"]:
            fmt = utils.format_openbabel(f"{fname}.{ext}")
            assert fmt == ext

    # .xyz extension has to be capitalized (XYZ)
    for fname in ["test", "root/test", "root/test.test"]:
        fmt = utils.format_openbabel(f"{fname}.xyz")
        assert fmt == "XYZ"

def test_deg_to_rad():

    assert utils.deg_to_rad(0) == pytest.approx(0)
    assert utils.deg_to_rad(90) == pytest.approx(np.pi / 2)
    assert utils.deg_to_rad(180) == pytest.approx(np.pi)
    assert utils.deg_to_rad(270) == pytest.approx(3 * np.pi / 2)
    assert utils.deg_to_rad(360) == pytest.approx(2 * np.pi)