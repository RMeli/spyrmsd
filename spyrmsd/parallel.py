try:
    from pebble import ProcessPool
except ImportError:
    errmsg = (
        "Parallel execution of SPyRMSD (`prmsdwrapper`) requires `pebble`."
        "Please install `pebble`: https://github.com/noxdafox/pebble"
    )
    raise ImportError(errmsg)

import os
from concurrent.futures import TimeoutError
from functools import partial
from typing import List, Optional, Union

import numpy as np

from spyrmsd import molecule
from spyrmsd.rmsd import rmsdwrapper


def prmsdwrapper(
    molrefs: Union[molecule.Molecule, List[molecule.Molecule]],
    mols: Union[molecule.Molecule, List[molecule.Molecule]],
    num_workers: Union[int, None] = 1,
    symmetry: bool = True,
    center: bool = False,
    minimize: bool = False,
    strip: bool = True,
    cache: bool = True,
    timeout: Optional[float] = None,
    chunksize: int = 1,
) -> List[float]:
    """
    Compute RMSD between two molecules with a timeout.
    Parameters
    ----------
    molrefs: Union[molecule.Molecule, List[molecule.Molecule]]
        Reference molecule
    mols: Union[molecule.Molecule, List[molecule.Molecule]]
        Molecules to compare to reference molecule
    num_workers: int
        Amount of processor to use for the parallel calculations
    symmetry: bool, optional
        Symmetry-corrected RMSD (using graph isomorphism)
    center: bool, optional
        Center molecules at origin
    minimize: bool, optional
        Minimised RMSD (using the quaternion polynomial method)
    strip: bool, optional
        Strip hydrogen atoms
    timeout: float, optional
        After how many seconds to stop the RMSD calculations
    chunksize: int, optional
        How many molecules to handle per process
    Returns
    -------
    List[float]
        RMSDs
    """

    # Ensure the num_workers is less or equal than the max number of CPUs.
    # Silencing MyPy since os.cpu_count() can return None
    num_workers = min(num_workers, os.cpu_count()) if os.cpu_count() is not None else 1  # type: ignore[type-var]

    # Cast the molecules to lists if they aren't already
    if not isinstance(molrefs, list):
        molrefs = [molrefs]
    if not isinstance(mols, list):
        mols = [mols]

    # Match the length of the molref
    if len(molrefs) == 1 and len(molrefs) < len(mols):
        molrefs = molrefs * len(mols)

    # Ensure molrefs and mols have the same len
    if not len(molrefs) == len(mols):
        raise ValueError("The input mol lists have different lengths")

    outputList = []

    timeoutCounter = 0
    errorCounter = 0

    with ProcessPool(max_workers=num_workers) as pool:
        rsmd_partial = partial(
            rmsdwrapper,
            symmetry=symmetry,
            center=center,
            minimize=minimize,
            strip=strip,
            cache=cache,
        )

        future = pool.map(
            rsmd_partial, molrefs, mols, timeout=timeout, chunksize=chunksize
        )
        iterator = future.result()
        idx = 0
        while True:
            try:
                result = next(iterator)
                outputList.append(result[0])
            except StopIteration:
                break
            except TimeoutError:
                # Not sure what the best way is to report the errors when the chunksize is non-zero
                timeoutCounter += 1  # * chunksize
                outputList.append(np.nan)
            except Exception as error:
                errorCounter += 1  # * chunksize
                outputList.append(np.nan)
                print(f"molecule pair {idx}: {error}")
            idx += 1

    print(
        f"{timeoutCounter} compounds timed out and were skipped, {errorCounter} compounds raised an error"
    )
    return outputList
