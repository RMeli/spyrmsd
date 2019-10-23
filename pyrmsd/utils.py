import numpy as np
import os


def format(fname: str) -> str:
    name, ext = os.path.splitext(fname)

    return ext[1:]  # Remove "."


def format_openbabel(fname: str) -> str:

    ext = format(fname)

    if ext == "xyz":
        ext = ext.upper()

    return ext


def deg_to_rad(angle: float) -> float:

    return angle * np.pi / 180.0


def rotate(v, angle, axis, units: str = "rad"):

    assert len(axis) == 3

    # Ensure rotation axis is normalised
    axis = axis / np.linalg.norm(axis)

    if units.lower() == "rad":
        pass
    elif units.lower() == "deg":
        angle = deg_to_rad(angle)
    else:
        raise Exception  # TODO

    t1 = axis * np.dot(axis, v)
    t2 = np.cos(angle) * np.cross(np.cross(axis, v), axis)
    t3 = np.sin(angle) * np.cross(axis, v)

    return t1 + t2 + t3
