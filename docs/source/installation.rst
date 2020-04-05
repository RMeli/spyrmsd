Installation
============

``spyrmsd`` is available on PyPI_ and conda-forge_ and can be easily installed from source.


Installation
------------

Installing ``spyrmsd`` with ``pip``, ``conda`` or from source will install the package as a library. In order to install the package as a standalone tool, `Open Babel`_ or RDKit_ need to be installed as well (see :ref:`Dependencies`).

pip
~~~

.. code-block:: bash

   pip install spyrmsd

conda
~~~~~

.. code-block:: bash

   conda install spyrmsd -c conda-forge

GitHub
~~~~~~

.. code-block:: bash

   git clone https://github.com/RMeli/spyrmsd.git
   cd spyrmsd
   pip install .

.. _Dependencies:

Dependencies
------------

``spyrmsd`` can be used both as a module or as a standalone tool.

Module
~~~~~~

The following packages are required to use ``spyrmsd`` as a module:

* graph-tool_ or NetworkX_
* numpy_
* scipy_
* qcelemental_

.. note::
   ``spyrmsd`` uses graph-tool_ by default but will  fall back  to NetworkX_ if the former is not installed (e.g. on Windows).

Standalone Tool
~~~~~~~~~~~~~~~

Additionally, one of the following packages is required to use ``spyrmsd`` as a standalone tool:

* `Open Babel`_
* RDKit_

.. note::
   RDKit_ has to be installed using ``conda`` (see `Why the RDKit isn't available on PyPi`_).

.. _PyPI: https://pypi.org/project/spyrmsd/
.. _conda-forge: https://github.com/conda-forge/spyrmsd-feedstock
.. _RDKit: https://rdkit.org/
.. _Why the RDKit isn't available on PyPi: https://rdkit.blogspot.com/2019/11/why-rdkit-isnt-available-on-pypi.html
.. _Open Babel: http://openbabel.org/
.. _graph-tool: https://graph-tool.skewed.de/
.. _NetworkX: https://networkx.github.io/
.. _numpy: https://numpy.org/
.. _scipy: https://www.scipy.org/
.. _qcelemental: http://docs.qcarchive.molssi.org/projects/qcelemental/en/latest/