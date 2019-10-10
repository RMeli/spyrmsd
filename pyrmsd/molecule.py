from pyrmsd import utils

from openbabel import pybel

def load(fname: str):

    fmt = utils.format_openbabel(fname)

    mol = next(pybel.readfile(fmt, fname))

    return  mol

