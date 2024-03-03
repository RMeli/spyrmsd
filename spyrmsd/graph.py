from spyrmsd import constants

import numpy as np
import warnings

_available_backends = []
_current_backend = None

## Backend aliases
_graph_tool_aliases = ["graph_tool", "graphtool", "graph-tool", "graph tool", "gt"]
_networkx_aliases = ["networkx", "nx"]

## Construct the alias dictionary
_alias_backendDict = {}
for alias in _graph_tool_aliases:
    _alias_backendDict[alias.lower()] = "graph-tool"
for alias in _networkx_aliases:
    _alias_backendDict[alias.lower()] = "networkx"

try:
    from spyrmsd.graphs.gt import (
        cycle as gt_cycle,
        graph_from_adjacency_matrix as gt_graph_from_adjacency_matrix,
        lattice as gt_lattice,
        match_graphs as gt_match_graphs,
        num_edges as gt_num_edges,
        num_vertices as gt_num_vertices,
        vertex_property as gt_vertex_property,
    )
    _available_backends.append("graph-tool")
except ImportError:
    warnings.warn("The graph-tool backend does not seem to be installed.", stacklevel=2)
    
try:
   from spyrmsd.graphs.nx import (
        cycle as nx_cycle,
        graph_from_adjacency_matrix as nx_graph_from_adjacency_matrix,
        lattice as nx_lattice,
        match_graphs as nx_match_graphs,
        num_edges as nx_num_edges,
        num_vertices as nx_num_vertices,
        vertex_property as nx_vertex_property,
    )
   _available_backends.append("networkx")
except ImportError:
   warnings.warn("The networkx backend does not seem to be installed.", stacklevel=2)


def _validate_backend(backend):
    standardized_backend = _alias_backendDict.get(backend.lower())
    if standardized_backend is None:
        raise ValueError("This backend is not recognized or supported")
    if standardized_backend not in _available_backends:
        raise ImportError(f"The {backend} backend doesn't seem to be installed")
    return standardized_backend

def available_backends():
    return _available_backends 

def set_backend(backend):
    ## Check if the backend is valid
    backend = _validate_backend(backend)

    ## Check if we actually need to switch backends
    if backend == _current_backend:
        warnings.warn(f"The backend is already {backend}", stacklevel=2)
        return _current_backend
        
    global cycle, graph_from_adjacency_matrix, lattice, match_graphs, num_edges, num_vertices, vertex_property
    
    if backend == "graph-tool":      
        cycle = gt_cycle
        graph_from_adjacency_matrix = gt_graph_from_adjacency_matrix
        lattice = gt_lattice
        match_graphs = gt_match_graphs
        num_edges = gt_num_edges
        num_vertices = gt_num_vertices
        vertex_property = gt_vertex_property
              
    elif backend == "networkx":
        cycle = nx_cycle
        graph_from_adjacency_matrix = nx_graph_from_adjacency_matrix
        lattice = nx_lattice
        match_graphs = nx_match_graphs
        num_edges = nx_num_edges
        num_vertices = nx_num_vertices
        vertex_property = nx_vertex_property      

    return backend

if len(_available_backends) == 0:
    raise ImportError("No valid backends found. Please ensure that either graph-tool or NetworkX are installed.")
else:
    if _current_backend is None:
        ## Set the backend to the first available (preferred) backend         
        _current_backend = set_backend(backend=_available_backends[0])
    
def get_backend():
    return _current_backend


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
