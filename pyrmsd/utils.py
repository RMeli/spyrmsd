import os

def format(fname: str) -> str:
    name, ext = os.path.splitext(fname)

    return ext[1:] # Remove "."

def format_openbabel(fname: str) -> str:

    ext = format(fname)

    if ext == "xyz":
        ext = ext.upper()

    return ext
