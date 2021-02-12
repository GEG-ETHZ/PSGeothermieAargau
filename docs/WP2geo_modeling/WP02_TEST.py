"""
Export a GemPy Model
======================

nananana NaNaNaNa EEEEYOOOO Gooodbye. ``connect`` of 'OpenWF'.
"""
# # Export a geological model from GemPy to use in MOOSE
# _implemented by [Jan Niederau](https://github.com/Japhiolite)_
# 
# This is a small example notebook guiding you through the process of exporting a geological model generated in [GemPy](https://www.gempy.org/) (Tutorial Chapter 1-1 therein) so it is usable as a Mesh in the [MOOSE](https://mooseframework.org/) framework.  
# 
# 
# These two lines are necessary only if GemPy is not installed 

import gempy as gp

import matplotlib.pyplot as plt

#%%
# Creating a geological model  
# ---------------------------
#
# The procedure of generating a geological model is presented in detail in [Chapter 1-1](https://nbviewer.jupyter.org/github/cgre-aachen/gempy/blob/master/notebooks/tutorials/ch1-1_Basics.ipynb) of the GemPy tutorials, so it will only be briefly presented here

geo_model = gp.create_model('tutorial_moose_exp')

gp.init_data(geo_model, [0,2000., 0,2000., 0,2000.], [50, 50, 80],
            path_o = "../../data/GemPy/simple_fault_model_orientations.csv",
            path_i = "../../data/GemPy/simple_fault_model_points.csv",
            default_values = True)

#%%
# present the units and series

geo_model.surfaces


#%%
# combine units in series and make two series, as the fault needs its own
gp.map_series_to_surfaces(geo_model,
                         {"Fault_Series" : 'Main_Fault',
                          "Strat_Series" : ('Sandstone_2', 'Siltstone', 'Shale', 'Sandstone_1', 'basement')},
                         remove_unused_series=True)

# set the fault series to be fault object
geo_model.set_is_fault(['Fault_Series'], change_color=False)

#%%
# check whether series were assigned correctly

geo_model.surfaces

#%%
# Model generation
# ----------------
# After loading in the data, we set it up for interpolation and compute the model.

gp.set_interpolator(geo_model,
                         compile_theano=True,
                         theano_optimizer='fast_compile',
                         verbose=[])

gp.compute_model(geo_model, compute_mesh=False)


gp.plot_2d(geo_model, direction='y', cell_number=45,show_data=False, show_boundaries=False, show_topography=False)


