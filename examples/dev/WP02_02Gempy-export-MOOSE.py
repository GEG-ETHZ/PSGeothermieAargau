#!/usr/bin/env python
# coding: utf-8

# # Export a geological model from GemPy to use in MOOSE
# _implemented by [Jan Niederau](https://github.com/Japhiolite)_
# 
# This is a small example notebook guiding you through the process of exporting a geological model generated in [GemPy](https://www.gempy.org/) (Tutorial Chapter 1-1 therein) so it is usable as a Mesh in the [MOOSE](https://mooseframework.org/) framework.  
# 
# 

# In[1]:


# These two lines are necessary only if GemPy is not installed 
import sys, os
sys.path.append("../../../gempy/")

import gempy as gp

import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# ## Creating a geological model  
# 
# The procedure of generating a geological model is presented in detail in [Chapter 1-1](https://nbviewer.jupyter.org/github/cgre-aachen/gempy/blob/master/notebooks/tutorials/ch1-1_Basics.ipynb) of the GemPy tutorials, so it will only be briefly presented here

# In[2]:


# Initiate a model
geo_model = gp.create_model('tutorial_moose_exp')


# In[3]:


gempy_dir = '../../../gempy/notebooks'
# Import data from CSV-files with setting the resolution and model extent
gp.init_data(geo_model, [0,2000., 0,2000., 0,2000.], [50, 50, 80],
            path_o = gempy_dir+"/data/input_data/tut_chapter1/simple_fault_model_orientations.csv",
            path_i = gempy_dir+"/data/input_data/tut_chapter1/simple_fault_model_points.csv",
            default_values = True);


# In[4]:


# present the units and series
geo_model.surfaces


# In[5]:


# combine units in series and make two series, as the fault needs its own
gp.map_series_to_surfaces(geo_model,
                         {"Fault_Series" : 'Main_Fault',
                          "Strat_Series" : ('Sandstone_2', 'Siltstone', 'Shale', 'Sandstone_1', 'basement')},
                         remove_unused_series=True);

# set the fault series to be fault object
geo_model.set_is_fault(['Fault_Series'], change_color=False)


# In[6]:


# check whether series were assigned correctly
geo_model.surfaces


# ## Model generation
# After loading in the data, we set it up for interpolation and compute the model.

# In[7]:


# set up interpolator
gp.set_interpolation_data(geo_model,
                          compile_theano=True, 
                          theano_optimizer='fast_compile',
                          verbose=[])


# In[8]:


# compute the model
gp.compute_model(geo_model, compute_mesh=False);


# In[9]:


# have a look at the data and computed model
gp.plot.plot_data(geo_model, direction='y')


# In[10]:


gp.plot.plot_section(geo_model, cell_number=24, direction='y',
                     show_data=False, show_legend=True)


# ## Exporting the Model to MOOSE  
# The voxel-model above already is the same as a model discretized in a hexahedral grid, so my immediately be used as input in a simulation tool, e.g. [MOOSE](https://mooseframework.org/). For this, we need to access to the unit IDs assigned to each voxel in GemPy. The array containing these IDs is called `lith_block`. 

# In[11]:


ids = geo_model.solutions.lith_block
print(ids)


# This array has the shape of `(x,)` and would be immediately useful, if GemPy and the chosen simulation code would _populate_ a grid in the same way. Of course, however, that is not the case. This is why we have to restructure the `lith_block` array, so it can be read correctly by MOOSE.

# In[12]:


# model resolution
nx, ny, nz = geo_model.grid.regular_grid.resolution

# model extent
xmin, xmax, ymin, ymax, zmin, zmax = geo_model.grid.regular_grid.extent


# These two parameters are important to, a) restructure `lith_block`, and b) write the input file for MOOSE correctly. 
# For a), we need to reshape `lith_block` again to its three dimensions and _re-flatten_ it in a _MOOSE-conform_ way.

# In[13]:


# reshape to 3D array
units = ids.reshape((nx, ny, nz))
# flatten MOOSE conform
units = units.flatten('F')


# The importance of `nx, ny, nz` is apparent from the cell above. But what about `xmin`, ..., `zmax`?  
# A MOOSE input-file for mesh generation has the following syntax:  
# 
# ```python
# [MeshGenerators]
#   [./gmg]
#     type = GeneratedMeshGenerator
#     dim = 3
#     nx = 50
#     ny = 50
#     nz = 80
#     xmin = 0.0
#     xmax = 2000.0
#     yim = 0.0
#     ymax = 2000.0
#     zmin = 0.0
#     zmax = 2000.0
#     block_id = '1 2 3 4 5 6'
#     block_name = 'Main_Fault Sandstone_2 Siltstone Shale Sandstone_1 basement'
#   [../]
# 
#   [./subdomains]
#     type = ElementSubdomainIDGenerator
#     input = gmg
#     subdomain_ids = ' ' # here you paste the transformed lith_block vector
#   [../]
# []
# 
# [Mesh]
#   type = MeshGeneratorMesh
# []
# ```
# 
# So these parameters are required inputs in the `[MeshGenerators]` object in the MOOSE input file. `GemPy` has a method to directly create such an input file, stored in `gempy.utils.export.py`.  
# 
# The following cell shows how to call the method:

# In[14]:


import gempy.utils.export as export
export.export_moose_input(geo_model, path='')


# This method automatically stores a file `geo_model_units_moose_input.i` at the specified path. Either this input file could be extended with parameters to directly run a simulation, or it is used just for creating a mesh. In the latter case, the next step would be, to run the compiled MOOSE executable witch the optional flag `--mesh-only`.  
# 
# E.g. with using the [PorousFlow module](https://mooseframework.inl.gov/modules/porous_flow/):
# 
# ```bash
# $path_to_moose/moose/modules/porous_flow/porous_flow-opt -i pct_voxel_mesh.i --mesh-only
# ```
# 
# How to compile MOOSE is described in their [documentation](https://mooseframework.inl.gov/getting_started/index.html). 
# 
# The now generated mesh with the name `geo_model_units_moose_input_in.e` can be used as input for another MOOSE input file, which contains the main simulation parameters. To call the file with the grid, the following part has to be added in the MOOSE simulation input file:  
# 
# ```python
# [Mesh]
#   file = geo_model_units_moose_input_in.e
# []
# ```
# 
# <hr>
# 
# The final output of the simulation may also be such an `.e`, which can, for instance, be opened with [paraview](https://www.paraview.org/). A simulated temperature field (purely conductive) of the created model would look like this:  
# 
# ![gempy_temperature](https://raw.githubusercontent.com/Japhiolite/a-Moose-and-you/master/imgs/GemPy_model_combined.png)
