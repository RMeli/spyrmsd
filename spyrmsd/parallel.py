try:
    from pebble import ProcessPool
except ImportError:
    errmsg = (
        "Parallel execution of SPyRMSD (`prmsdwrapper`) requires `pebble`."
        "Please install `pebble`: https://github.com/noxdafox/pebble"
    )
    raise ImportError(errmsg)

import os
import warnings
from concurrent.futures import TimeoutError
from functools import partial
from typing import List, Optional, Union

import numpy as np

from spyrmsd import molecule
from spyrmsd.rmsd import rmsdwrapper


def prmsdwrapper(
    molrefs: Union[molecule.Molecule, List[molecule.Molecule]],
    mols: Union[molecule.Molecule, List[molecule.Molecule]],
    symmetry: bool = True,
    center: bool = False,
    minimize: bool = False,
    strip: bool = True,
    cache: bool = True,
    num_workers: Optional[int] = None,
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
    if num_workers is None:
        num_workers = os.cpu_count()
    num_workers = min(num_workers, os.cpu_count())  # type: ignore[type-var]

    if chunksize > 1 and timeout is not None:
        warnings.warn("When using the timeout feature, a chunksize of 1 is required")
        ## We could also raise an error here
        chunksize = 1

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
        raise ValueError("The 'mols' and 'molrefs' lists have different lengths.")

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
                timeoutCounter += 1
                outputList += [np.nan] * chunksize
            except Exception:
                errorCounter += 1
                outputList.append(np.nan)
            idx += 1
    if timeoutCounter + errorCounter > 0:
        # Calculate total number of np.nan
        failedCompoundsTotal = np.count_nonzero(np.isnan(outputList))

        warnings.warn(
            f"{failedCompoundsTotal} compounds failed to process successfully ({errorCounter} compounds raised an error, {timeoutCounter} chunks timed out"
        )
    return outputList
