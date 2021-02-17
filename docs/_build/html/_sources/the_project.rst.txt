The project
===========

The project **Pilot Study Geothermie Aargau** is a three year project conducted at the Geothermal Energy and Geofluids Group (ETH Zurich), in cooperation with 
swisstopo and the Georessources Switzerland Group (also ETH). Four main work packages constitute to the pilot study, ranging from data acquisition and analysis, over
geological modeling to geothermal reservoir simulation. 

All these steps should be fully reproducible, and thus are used with open-source software. To be fully reproducible, the project relies on an extensive documentation, which
this site is part of. Presented here are main examples from each of the main work packages to describe the process of getting from diverse data sources to a model-based heat flow estimate at 
variable depths. 

The project consists of 4 main work packages, which are shown in the figure below. Each work package is connected to one main software or language, used within this work package. 

-  Work package one comprises data collection and compilation in SQL databases (using sqlite)   
-  Work package two focuses on geological modeling of a geothermal reservoir system, i.e. getting the overall structure of the reservoir system   
-  Work package three then uses the results of work package two for getting to heat-transport simulations using open source codes, such as MOOSE or SHEMAT-Suite   
-  Work package four comprises package and workflow documentation. As such, this docs-page is part of WP4, which heavily relies on Sphinx and Jupyter Notebooks.  

.. image:: ./_static/logos/ps_aargau_project_WPs.png
  :width: 800
  :alt: Graphical description of the 4 work packages
  :class: with-shadow 


.. autoclass:: OpenWF.db_access.connect

.. include:: gen_modules/backreferences/OpenWF.db_access.connect.examples
.. raw:: html

    <div class="sphx-glr-clear"></div>
