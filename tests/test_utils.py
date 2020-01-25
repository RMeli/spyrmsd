import random

import numpy as np
import pytest

from spyrmsd import utils


@pytest.mark.parametrize("ext", ["pdf", "smi", "xyz"])
def test_format(ext: str) -> None:
    for fname in ["test", "root/test", "root/test.test"]:
        fmt = utils.format(f"{fname}.{ext}")
        assert fmt == ext


@pytest.mark.parametrize(
    "extin, extout", [("smi", "smi"), ("pdb", "pdb"), ("xyz", "XYZ")]
)
def test_molformat(extin: str, extout: str) -> None:

    for fname in ["test", "root/test", "root/test.test"]:
        fmt = utils.molformat(f"{fname}.{extin}")
        assert fmt == extout


@pytest.mark.parametrize(
    "deg, rad",
    [(0, 0), (90, np.pi / 2), (180, np.pi), (270, 3 * np.pi / 2), (360, 2 * np.pi)],
)
def test_deg_to_rad(deg: float, rad: float) -> None:

    assert utils.deg_to_rad(deg) == pytest.approx(rad)


def test_rotate_invalid():

    with pytest.raises(ValueError):
        utils.rotate(np.array([1, 0, 0]), 0, np.array([0, 0, 1]), units="none")


@pytest.mark.parametrize("z", [np.array([0, 0, random.random()]) for _ in range(10)])
def test_rotate_z(z: np.ndarray) -> None:

    for deg, rad in [(0, 0), (45, np.pi / 4), (90, np.pi / 2)]:
        v_deg = utils.rotate(np.array([1, 0, 0]), deg, z, units="deg")
        v_rad = utils.rotate(np.array([1, 0, 0]), rad, z, units="rad")

        assert np.allclose(v_deg, v_rad)
        assert np.allclose(v_deg, np.array([np.cos(rad), np.sin(rad), 0]))
