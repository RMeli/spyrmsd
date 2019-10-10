import numpy as np
import os

def format(fname: str) -> str:
    name, ext = os.path.splitext(fname)

    return ext[1:] # Remove "."

def format_openbabel(fname: str) -> str:

    ext = format(fname)

    if ext == "xyz":
        ext = ext.upper()

    return ext

def deg_to_rad(angle):

    return angle * np.pi / 180.0
