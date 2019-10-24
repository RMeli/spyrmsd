import numpy as np
import os


def format(fname: str) -> str:
    """
    Extract format extension from file name.

    Parameters
    ----------
    fname : str
        File name

    Returns
    -------
    str
        File extension

    Notes
    -----
    The file extension is returned without the `.` character, i.e. for the file
    `path/filename.ext` the string `ext` is returned.
    """
    name, ext = os.path.splitext(fname)

    return ext[1:]  # Remove "."


def format_openbabel(fname: str) -> str:
    """
    Extract an OpenBabel-friendly format from file name.

    Parameters
    ----------
    fname : str
        File name

    Returns
    -------
    str
        File extension in an OpenBabel-friendly format

    Notes
    -----
    File types in OpenBabel do not always correspond to the file extension. This
    function converts the file extension to an OpenBabel file type.

    The following table shows the different conversions performed by this function:

    ========= =========
    Extension File Type
    --------- ---------
    xyz       XYZ
    ========= =========
    """

    ext = format(fname)

    if ext == "xyz":
        # xyz files in OpenBabel are called XYZ
        ext = ext.upper()

    return ext


def deg_to_rad(angle: float) -> float:
    """
    Convert angle in degrees to angle in radians.

    Parameters
    ----------
    angle : float
        Angle (in degrees)

    Returns
    -------
    float
        Angle (in radians)
    """

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
