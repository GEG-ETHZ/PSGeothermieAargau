Getting started
===============

The code project ``OpenWF`` is still in active development. Current procedure for installation is to clone this repository
and import it manually to your projects. A more convenient installation via pypi may come in the future.

Cloning the repository and adding it to the path, i.e. installation
-------------------------------------------------------------------

In a first step, clone this repository using the following command (or by manually downloading the zip file from the GitHub page)
(via HTTPS)

.. code-block:: bash

   $ git clone https://github.com/GEG-ETHZ/PSGeothermieAargau.git

(via SSH)
.. code-block:: bash

   $ git clone git@github.com:GEG-ETHZ/PSGeothermieAargau.git


using the python package ``sys`` you then append the path to the repository

.. code-block:: python    
   
   import sys
   sys.path.append("<path/to/cloned/repository/PSGeothermieAargau>")


followed by the import of the module

.. code-block:: python

   import OpenWF as owf


This package has been tested with Python 3.8 and works best currently with this pyhton Version.

``OpenWF`` depends on a couple of other Python libraries:

* ``GemPy`` for creating the geological models
* ``PyMC3`` for bayesian inference (_optional_) Without PyMC3, only Massive Monte Carlo is possible.
* ``matplotlib`` for plotting
* ``numpy`` for efficient numerical methods
* ``sqlite`` for database access and manipulation  
* ...
  
an ``environment.yml`` file is provided for installation on Windows and Linux. These files can be used for creating
functioning anaconda environments for working with ``OpenWF``.


Using the package
-----------------

Within a python script or a jupyter notebook, append the path to the cloned repository before using the methods. As an
example, we show how one would calculate the conductive heat flow from an HDF5 file:

.. code-block:: python
   import sys
   sys.path.append('path/to/PSGeothermieAargau/repository/')
   import OpenWF as owf
   import OpenWF.postprocessing as opp

   model_path = 'path/to/model/outputs/in/HDF5/format.h5'
   fid = pp.read_hdf_file(model_path, write=False)
   qx, qy, qz = opp.calc_cond_hf(fid)
   print('calc_cond_hf yields conductive heat flow in all three model dimensions')
