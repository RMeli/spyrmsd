from pyrmsd import utils

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