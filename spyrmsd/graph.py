import numpy as np
import sys

from spyrmsd import constants

## Checking and initializing the available backends
graph_backends = []

try:
    import spyrmsd.graphs.gt
    graph_backends.append("graph_tool")

except ImportError:
    print("Graph Tool backend not found")
    
try:
   import spyrmsd.graphs.nx
   graph_backends.append("networkx")
   
except ImportError:
    print("NetworkX backend not found")

if len(graph_backends) == 0:
    sys.exit("No valid graph backends were found, please make sure atleast one of the supported backends is installed correctly")

def get_backends():
    return graph_backends

def adjacency_matrix_from_atomic_coordinates(
    aprops: np.ndarray, coordinates: np.ndarray
) -> np.ndarray:
    """
    Compute adjacency matrix from atomic coordinates.

    Parameters
    ----------
    aprops: numpy.ndarray
        Atomic properties
    coordinates: numpy.ndarray
        Atomic coordinates

    Returns
    -------
    numpy.ndarray
        Adjacency matrix

    Notes
    -----

    This function is based on an automatic bond perception algorithm: two
    atoms are considered to be bonded when their distance is smaller than
    the sum of their covalent radii plus a tolerance value. [3]_

    .. warning::
        The automatic bond perceptron rule implemented in this functions
        is very simple and only depends on atomic coordinates. Use
        with care!

    .. [3] E. C. Meng and R. A. Lewis, *Determination of molecular topology and atomic
       hybridization states from heavy atom coordinates*, J. Comp. Chem. **12**, 891-898
       (1991).
    """

    n = len(aprops)

    assert coordinates.shape == (n, 3)

    A = np.zeros((n, n))

    for i in range(n):
        r_i = constants.anum_to_covalentradius[aprops[i]]

        for j in range(i + 1, n):
            r_j = constants.anum_to_covalentradius[aprops[j]]

            distance = np.sqrt(np.sum((coordinates[i] - coordinates[j]) ** 2))

            if distance < (r_i + r_j + constants.connectivity_tolerance):
                A[i, j] = A[j, i] = 1

    return A
