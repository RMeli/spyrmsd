xyzgraph Adapter
================

``spyrmsd`` can use `xyzgraph <https://github.com/aligfellow/xyzgraph>`_ as an
optional adjacency matrix builder for molecular graph construction.

Install the optional dependency with:

.. code-block:: bash

   pip install spyrmsd[xyzgraph]

The adapter keeps the same high-level workflow already used by ``spyrmsd``:
first build an adjacency matrix, then pass it to :func:`spyrmsd.rmsd.symmrmsd`.

Example
-------

.. code-block:: python

   from spyrmsd.adapters.xyzgraph import (
       adjacency_matrix_from_atomic_coordinates as xyzgraph_adjacency_matrix,
   )
   from spyrmsd.rmsd import symmrmsd

   amref = xyzgraph_adjacency_matrix(atomicnums_ref, coords_ref)
   am = xyzgraph_adjacency_matrix(atomicnums, coords)

   value = symmrmsd(
       coords_ref,
       coords,
       atomicnums_ref,
       atomicnums,
       amref,
       am,
   )

.. note::
   The default graph builder remains the built-in ``simple`` one,
   based on distance and VdW radii. The ``xyzgraph`` adapter is opt-in.
