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


def test_rotate_invalid():

    with pytest.raises(Exception):
        utils.rotate(np.array([1, 0, 0]), 0, np.array([0, 0, 1]), units="none")


def test_rotate_z():

    for z in [np.array([0, 0, 1]), np.array([0, 0, 2])]:

        for deg, rad in [(0, 0), (45, np.pi / 4), (90, np.pi / 2)]:
            v_deg = utils.rotate(np.array([1, 0, 0]), deg, z, units="deg")
            v_rad = utils.rotate(np.array([1, 0, 0]), rad, z, units="rad")

            assert np.allclose(v_deg, v_rad)
            assert np.allclose(v_deg, np.array([np.cos(rad), np.sin(rad), 0]))
