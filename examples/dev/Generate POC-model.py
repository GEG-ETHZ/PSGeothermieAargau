#!/usr/bin/env python
# coding: utf-8

# # Create a 3D model of a Permo-Carboniferous Trough (PCT)
# 
# Based on a seismic section from the NAGRA report [NAGRA NAB 14-17](https://www.nagra.ch/data/documents/database/dokumente/$default/Default%20Folder/Publikationen/NABs%202004%20-%202015/d_nab14-017.pdf), we extracted interface and orientation points for lithological units and faults.  
# 
# The lithological units comprise the permo-carboniferous filling (divided in three stages based on the report results), Mesozoic, Tertiary, and Quaternary strata, as well as the Palaeozoic crystalline basement rocks.

# In[1]:


import warnings
warnings.filterwarnings("ignore")

# Importing GemPy
import gempy as gp
from gempy.plot import visualization_2d as vv

# Importing auxilary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# For embedding matplotlib figures
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


print(f"Code run with GemPy version: {gp.__version__}")


# # Initialize the model
# We start with modelling the trough by generating a gempy model object. This will use interface points and orientations, which we previously stored in a `.csv` file.

# In[3]:


np.random.seed(333)
# Import data
# Create a model instance
geo_model = gp.create_model('POC_model')

# Initialize the model, set dimension and load interface and orientation data
gp.init_data(geo_model, [0, 28000., 0, 14000., -6500, 1000.], [100, 50, 60],
            path_i = '../data/leu2008/line82_interfaces_wo_middle_MC.csv',
            path_o = '../data/leu2008/line82_foliations_wo_middle_MC.csv');
geo_model.set_topography(source='random', d_z=np.array([300,1000]));


# In[4]:


gp.plot_2d(geo_model, show_data=True, show_topography=True);


# In[5]:


nx, ny, nz = geo_model.grid.regular_grid.resolution
xmin, xmax, ymin, ymax, zmin, zmax = geo_model.solutions.grid.regular_grid.extent

delx = xmax / nx
dely = ymax / ny


# In[6]:


bhole = np.array([[32, 15],
                 [79, 23],
                 [53, 34]])


# In[7]:


gp.plot_2d(geo_model, show_data=False, show_topography=True, show_values=False, section_names=['topography'],
          show_results=False, show_legend=False)
plt.scatter(bhole[:,0]*delx, bhole[:,1]*dely, s=500, marker='^', facecolor='k')
plt.show()


# # Adding information to the model
# ## Surfaces

# In[4]:


col_dict = {'basement': '#c7848f',
           'Lower-filling': '#a5d490', 
           'Upper-filling': '#cfc199',
           'Unconformity': '#725c9a',
           'Orange': '#ff792b',
           'Pink': '#e588f3',
           'Tertiary': '#dbdbac',
           'Fault2': '#015482',
           'Fault5': '#015482',
           'Fault6': '#015482',
           'Thrust1_south': '#5DA629',
           'Thrust2_south': '#5DA629'}
geo_model.surfaces.colors.change_colors(col_dict)
geo_model.surfaces


# ## Characteristics  
# Main features of the model is the asymetric graben system, with the major fault (denoted with **A**), and the graben fill, which is not present beyond the graben shoulders. This, as well as the stop of major faults beneath the mesozoic units (blue units) are important considerations for the modelling process.  
# These could be caught, for instance, in likelihood functions if we model the PCT as a Bayesian inference problem.

# In[5]:


# Assign formations to series
gp.map_series_to_surfaces(geo_model,
                         {"Thrust1_series": 'Thrust1_south',
                          "Thrust2_series": 'Thrust2_south',
                          "Fault2_series": 'Fault2',
                          "Fault5_series": 'Fault5',
                          "Fault6_series": 'Fault6',
                         "Post_tectonic_series": ('Tertiary', 'Pink', 'Orange'),
                          "Detachement": 'Unconformity',
                         "Syn_tectonic_series2": 'Upper-filling',
                         #"Syn_tectonic_series1": 'Middle-filling',
                         "Pre_tectonic_series": 'Lower-filling'},
                         remove_unused_series=True)
geo_model.surfaces


# In[6]:


# Set Faults
geo_model.set_is_fault(['Thrust1_series', 'Thrust2_series',
                        'Fault2_series', 'Fault5_series', 'Fault6_series'],
                      change_color=False);
#geo_model.set_is_finite_fault(series_fault=['Fault1_series', 'Fault7_series', 'Fault6_series', 
#                                            'Fault5_series', 'Fault3_series', 'Fault4_series'],
#                              toggle=True)


# In[7]:


geo_model.set_bottom_relation(series=['Post_tectonic_series', 
                                      'Pre_tectonic_series',
                                      #'Syn_tectonic_series1',
                                      'Syn_tectonic_series2'], bottom_relation='Onlap') #,


# In[8]:


# table of fault relations
geo_model.faults.faults_relations_df


# In[9]:


fr = np.array([[False, True, False, False, False, True, False, False,   False, False],
               [False, False, False, False, False, True, False, False,  False, False],
               [False, False, False, False, False, False, True, True,  True, True],
               [False, False, False, False, False, False, True, True,  True, True],
               [False, False, False, False, False, False, True, True,  True, True],
               [False, False, False, False, False, False, False, False, False, False],
               [False, False, False, False, False, False, False, False, False, False],
               [False, False, False, False, False, False, False, False, False, False],
               [False, False, False, False, False, False, False, False, False, False],
               [False, False, False, False, False, False, False, False, False, False]])
geo_model.set_fault_relation(fr)


# ## Creating the model
# Now that we set the parameters and fault relations, it is time to start the modeling process:

# In[10]:


# decrease the kriging range
geo_model.modify_kriging_parameters('range', 20000.)
geo_model.modify_kriging_parameters('$C_o$', 2e5)
#geo_model.modify_surface_points('all', smooth=1e-6)

# Set the interpolator function
# Create the theano model
gp.set_interpolator(geo_model,
                         compile_theano=True,
                         theano_optimizer='fast_compile',
                         verbose=[],
                         update_kriging=False);


# In[11]:


# Compute the model
sol = gp.compute_model(geo_model)


# In[12]:


geo_model.save_model(name='POC_PCT_model', 
                     path='../models/2021-06-04_POC_base_model')


# When plotting the scalar field, we actually see what the gradients (the orientation triangles) look like for the different units / faults. 

# In[18]:


gp.plot_2d(geo_model, cell_number=25, direction='y', show_data=False, show_topography=False,
          show_lith=True, show_results=True, show_boundaries=True)


# In[15]:


gp.plot_2d(geo_model, cell_number=50, direction='z', show_data=False, show_topography=False,
          show_lith=True, show_results=True, show_boundaries=False)


# ## Simulate Gravity
# Using the now generated POC-model, we simulate its gravity at different locations. These locations will be treated as observations later on in the workflow. 
# In a first step, we distribute 15 points randomly across the topography of our model. Those will be the station locations, where we pick up the gravity signal of the POC-model.

# In[28]:


# distribute stations
import random
np.random.seed(58)
station_indices = np.random.randint(0, high=4999, size=15)
station_coordinates = geo_model._grid.topography.values[station_indices, :]

cs = plt.scatter(station_coordinates[:,0], station_coordinates[:,1], c=station_coordinates[:,2], cmap='viridis')
plt.colorbar(cs)


# In[29]:


from gempy.assets.geophysics import GravityPreprocessing
geo_model.set_centered_grid(station_coordinates,  resolution = [10, 10, 15], radius=6000)
g = GravityPreprocessing(geo_model.grid.centered_grid)
tz = g.set_tz_kernel()


# In[30]:


# add densities - from abdelfettah 2014 and SAPHYR
densities = [0, 0, 0, 0, 0, 2.466, 2.61, 2.53, 
             2.61, 2.47, 2.55, 2.67]
geo_model.add_surface_values(densities, ['density'])


# In[31]:


geo_model.set_active_grid('centered', reset=True)


# In[32]:


gp.set_interpolator(geo_model, output=['gravity'], theano_optimizer='fast_run', update_kriging=False)
sol = gp.compute_model(geo_model)
grav = sol.fw_gravity


# In[24]:


grav1 = grav.reshape(len(grav),1)
station_forw_grav = np.append(station_coordinates, grav1, axis=1)


# In[25]:


np.savetxt('../models/20210322_forw_grav_seed58.csv', station_forw_grav, fmt='%.2f, %.2f, %.2f, %.5f')

