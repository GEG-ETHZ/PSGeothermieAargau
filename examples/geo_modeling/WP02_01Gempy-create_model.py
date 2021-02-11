#!/usr/bin/env python
# coding: utf-8

# # Create a geological model in GemPy
# Geological modeling is the core task of Work Package 2 in this project. We use the opensource code GemPy for generating geological models whose purpose is to provide structure for later heat-transport simulations. 
# 
# In this notebook, we demonstrate how to build a geological model using input data from the 2D seismic line 82 in Leu 2008 \cite{leu_permokarbon-kartenskizze_2008}

# In[1]:


import sys, os
sys.path.append("../..")

import gempy as gp
get_ipython().run_line_magic('matplotlib', 'inline')

import numpy as np
import pandas as pd


# In[2]:


# Import data
data_path = '../../data/processed/GemPy/'
geo_model = gp.create_model('Permo_Carb_Trough')
gp.init_data(geo_model, [0, 30000., 0, 10., -7000, 1000.], [100, 1, 100],
                         path_o = data_path+'line82_foliations.csv',
                         path_i = data_path+'line82_interfaces.csv', default_values=False);


# In[3]:


gp.get_data(geo_model, 'surface_points').head()


# In[4]:


import matplotlib.pyplot as plt
gp.plot.plot_data(geo_model, direction='y')


# In[5]:


import matplotlib.image as mimg
img = mimg.imread('../../figs/line82-nf-10.png')
plt.imshow(img);


# In[5]:


# Assign formations to series
gp.map_series_to_surfaces(geo_model,
                         {"Fault7_series": 'Fault7',
                          "Fault1_series": 'Fault1',
                          "Fault6_series": 'Fault6',
                          "Fault2_series": 'Fault2',
                          "Fault5_series": 'Fault5',
                         "Fault3_series": 'Fault3',
                         "Fault4_series": 'Fault4',
                         "Post_tectonic_series": ('Quaternary', 'Tertiary', 'Mesozoic'),
                         "Syn_tectonic_series": 'Upper-filling',
                         "Pre_tectonic_series": ('Lower-filling','Middle-filling')},
                         remove_unused_series=True)
geo_model.surfaces


# In[16]:


geo_model.set_is_fault(['Fault7_series', 'Fault1_series', 'Fault6_series', 'Fault2_series', 
                        'Fault5_series', 'Fault3_series', 'Fault4_series'])
#geo_model.set_is_finite_fault(series_fault=['Fault6_series', 'Fault5_series', 'Fault3_series', 'Fault4_series'],
#                              toggle=True)


# In[6]:


geo_model.faults.faults_relations_df


# In[7]:


gp.activate_interactive_df(geo_model)


# In[8]:


geo_model.qi.get('faults_relations')


# In[14]:


geo_model.faults.faults_relations_df


# In[15]:


gp.set_interpolation_data(geo_model,
                         compile_theano=True,
                         theano_optimizer='fast_compile',
                         verbose=[])
# OLD GEMPY
# interp_data = gp.InterpolatorData(geo_data,
#                                 output='geology', compile_theano=True,
#                                 theano_optimizer='fast_compile')


# In[16]:


sol = gp.compute_model(geo_model, compute_mesh=False)
# OLD GEMPY
#lith_block, fault_block = gp.compute_model(interp_data)


# In[17]:


gp.plot.plot_scalar_field(geo_model, cell_number=0, series=7)


# In[22]:


get_ipython().run_line_magic('matplotlib', 'inline')
gp.plot.plot_section(geo_model, cell_number=0, direction='y', show_faults=True,
                     show_data=True, interpolation='nearest', ve=1)


# In[8]:


gp.plotting.plot_section(geo_data, lith_block[0], cell_number=24, direction='y', plot_data=True)


# In[59]:


import matplotlib.image as mimg
img = mimg.imread('Screenshot 2019-02-13 at 15.14.59.png')
plt.imshow(img)


# In[15]:


geo_model.save_model(name='section_test_2D', path='../')


# In[16]:


ids = geo_model.solutions.lith_block
ids = ids.reshape((ids.shape[0],1)).astype(int)


# In[17]:


ids.shape


# In[18]:


grid = geo_model.solutions.grid.values


# In[19]:


xyz = np.hstack((grid,ids))


# In[20]:


np.savetxt('../section_test_2D/model_grid_2D', xyz, fmt="%.3f, %.3f, %.3f, %i",
          delimiter=",", header="#x,y,z,id", comments='#')


# In[34]:


grid[:-1:100,0]


# In[66]:


x, z = np.meshgrid(grid[:-1:100,0], grid[:40,2])


# In[67]:


ids_resh = ids.reshape(geo_model.grid.regular_grid.resolution[0],
                      geo_model.grid.regular_grid.resolution[1],
                      geo_model.grid.regular_grid.resolution[2]).T


# In[70]:


plt.pcolor(x, z, ids_resh[:,5,:], cmap='viridis')


# In[32]:


_x, _y, _z = (slice(0, geo_model.grid.regular_grid.resolution[0]),
             slice(0, geo_model.grid.regular_grid.resolution[1]),
             slice(0, geo_model.grid.regular_grid.resolution[2]))
extent_v = geo_model.grid.extent[[0, 1, 4, 5]]


# In[35]:


ids_resh[_x, _y, _z].T.shape


# In[27]:


extent_v


# In[71]:


np.savetxt('../section_test/lith_vector', ids, fmt="%i", header="id", comments="#")


# # References
# 
# (<a id="cit-leu_permokarbon-kartenskizze_2008" href="#call-leu_permokarbon-kartenskizze_2008">Leu, 2008</a>) Werner Leu, ``_Permokarbon-Kartenskizze (Rohstoffe). Kompilation eines GIS-Datensatzes auf der Basis von bestehenden Unterlagen (Bereich Schweizer Mittelland)_'',  2008.
# 
# 
