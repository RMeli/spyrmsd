Tutorial
========

.. code:: ipython3

    import spyrmsd
    from spyrmsd import io, rmsd

OpenBabel or RDKit
~~~~~~~~~~~~~~~~~~

``spyrmsd`` natively supports
`OpenBabel <http://openbabel.org/wiki/Main_Page>`__ and
`RDKit <http://rdkit.org/>`__ to load molecules in order to work as a
standalone tool. However, the API for RMSD calculations is extremely
minimal and only needs the following information:

-  Atomic coordinates
-  Atomic numbers
-  Molecular adjacency matrix (for symmetry)

This means that ``spyrmsd`` can be used in combination with any library
that can provide such information.

Loading Molecules
~~~~~~~~~~~~~~~~~

The ``spyrmsd.io`` module provides functions to easily load molecules
from a file and to transform them into a ``spyrmsd.molecule.Molecule``
object:

.. code:: ipython3

    ref = io.loadmol("molecules/1a4k_ligand.sdf")

``io.loadmol`` load a single molecule from a file. In order to load all
molecules we need to use ``io.loadallmols``:

.. code:: ipython3

    mols = io.loadallmols("molecules/1a4k_dock.sdf")

Loading RDKit or OpenBabel Molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``spyrmsd`` natively supports Open Babel and RDKit (if installed). The
``Molecule`` class provides ``from_openbabel()`` and ``from_rdkit()``
constructors.

.. code:: ipython3

    from rdkit import Chem
    from rdkit.Chem import AllChem

    rdmol1 = Chem.MolFromSmiles("c1ccccc1")
    rdmol2 = Chem.MolFromSmiles("c1ccccc1")
    AllChem.EmbedMolecule(rdmol1)
    AllChem.EmbedMolecule(rdmol2)

    from spyrmsd.molecule import Molecule
    from spyrmsd.rmsd import rmsdwrapper

    mol1 = Molecule.from_rdkit(rdmol1)
    mol2 = Molecule.from_rdkit(rdmol2)

    rmsdwrapper(mol1, mol2)


.. parsed-literal::

    <frozen importlib._bootstrap>:241: RuntimeWarning: to-Python converter for std::__1::pair<double, double> already registered; second conversion method ignored.
    [21:58:01] Molecule does not have explicit Hs. Consider calling AddHs()
    [21:58:01] Molecule does not have explicit Hs. Consider calling AddHs()




.. parsed-literal::

    [0.019162902039384797]



Removing Hydrogen Atoms
~~~~~~~~~~~~~~~~~~~~~~~

Hydrogen atoms can be removed with the ``strip()`` function:

.. code:: ipython3

    ref.strip()

.. code:: ipython3

    for mol in mols:
        mol.strip()

Symmetry-Corrected RMSD
~~~~~~~~~~~~~~~~~~~~~~~

``spyrmsd`` only needs atomic coordinates, atomic number and the
molecular adjacency matrix to compute the standard RMSD with
``spyrmsd.rmsd.symmrmsd``. The ``spyrmsd.molecule.Molecule`` class
provides easy access to such information:

.. code:: ipython3

    coords_ref = ref.coordinates
    anum_ref = ref.atomicnums
    adj_ref = ref.adjacency_matrix

.. code:: ipython3

    coords = [mol.coordinates for mol in mols]
    anum = mols[0].atomicnums
    adj = mols[0].adjacency_matrix

With this information we can easily compute the RMSD between the
reference molecule and all other molecules:

.. code:: ipython3

    RMSD = rmsd.symmrmsd(
        coords_ref,
        coords,
        anum_ref,
        anum,
        adj_ref,
        adj,
    )

    print(RMSD)


.. parsed-literal::

    [2.0246085732404446, 1.4951562971486378, 10.028009301306854, 7.900570020309068, 7.578344354783399, 9.52999506817054, 4.952371789159667, 7.762808670066815, 9.996922964463582, 7.1732072690335755]


Minimum RMSD
------------

We can also compute the minimum RMSD obtained by superimposing the
molecular structures:

.. code:: ipython3

    RMSD = rmsd.symmrmsd(
        coords_ref,
        coords,
        anum_ref,
        anum,
        adj_ref,
        adj,
        minimize=True,
    )

    print(RMSD)


.. parsed-literal::

    [1.2012368667355435, 1.0533413220699535, 1.153253104575529, 1.036542688936588, 0.8407673221224187, 1.1758143217869736, 0.7817315189656655, 1.0933314311267845, 1.0260767175206462, 0.9586369647000478]



Change Backend
~~~~~~~~~~~~~~

``spyrmsd`` supports multiple backends. You see which backends are
available by looking at the ``available_backends`` attribute:

.. code:: ipython3

    spyrmsd.available_backends




.. parsed-literal::

    ['graph_tool', 'networkx']



The available backends are a subset of the supported backends. Only the
backends that are installed will be available.

You can check the current backend with

.. code:: ipython3

    spyrmsd.get_backend()




.. parsed-literal::

    'graph_tool'



You can switch the backend using

.. code:: ipython3

    spyrmsd.set_backend("networkx")
    spyrmsd.get_backend()




.. parsed-literal::

    'networkx'
