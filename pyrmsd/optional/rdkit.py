import rdkit.Chem as Chem

from pyrmsd import utils


def load(fname: str) -> Chem.rdchem.Mol:
    """
    Load molecule from file

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    Chem.rdchem.Mol
        RDkit molecule
    """

    fmt = utils.molformat(fname)

    if fmt == "mol2":
        rdmol = Chem.MolFromMol2File(fname)

    return rdmol
