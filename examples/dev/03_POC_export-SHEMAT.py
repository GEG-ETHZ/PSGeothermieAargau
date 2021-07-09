"""
Create SHEMAT-Suite models
==========================
 
 With the MC ensemble of generated geological models stored in the respective lith-blocks, we can use them to create SHEMAT-Suite models for then doing 
 forward simulations of Heat- and Mass-transfer.
"""


# Libraries
import os,sys
import glob
import numpy as np
import itertools as it
import gempy as gp
import pandas as pd
print(f"Run with GemPy version {gp.__version__}")


# load the base model
# Gempy 2.2.9
model_path = '../../models/2021-04-01_GeoMol_Boreholes_Erode/PCT_erode_thick_fHorizonfloat64/'
geo_model = gp.load_model('PCT_erode_thick_fHorizon',
                         path=model_path, recompile=False)


# In[3]:


geo_model._grid.regular_grid.extent


# In[4]:


topo = geo_model._grid.topography.values.shape
topo_mask = geo_model._grid.regular_grid.mask_topo


# In[5]:


def topomask(geo_model, lith_block):
    """
    Method to mask the GemPy model with topography and change masked IDs to a new one for air.
    """
    model_shape = geo_model._grid.regular_grid.mask_topo.shape
    topo_mask = geo_model._grid.regular_grid.mask_topo
    
    # round lithblock and change to integers
    ids = np.round(lith_block)
    ids = ids.astype(int)
    
    mask_lith = ids.reshape(model_shape)
    mask_lith[topo_mask] = np.unique(mask_lith)[-1] + 1
    
    return mask_lith

def conc_lithblocks(path):
    """
    Load lithblocks npy files and concatenate them
    """
    fids = glob.glob(path+'*.npy')
    lith_blocks = []
    
    for f in fids:
        blocks = np.load(f)
        lith_blocks.append(blocks)
        all_lith_blocks = np.concatenate(lith_blocks)
        
    return all_lith_blocks


# In[6]:


def export_shemat_suite_input_file(geo_model, lithology_block, 
                                   units: pd.DataFrame=None, bcs_file: str=None, 
                                   data_file: str=None, borehole_logs: np.array=None,
                                   path: str=None, filename: str='geo_model_SHEMAT_input_erode'):
    """
    Method to export a 3D geological model as SHEMAT-Suite input-file for a conductive HT-simulation. 

    Args:
        path (str): Filepath for the exported input file (default './')
        filename (str): name of exported input file (default 'geo_model_SHEMAT_input')
    """
    # get model dimensions
    nx, ny, nz = geo_model.grid.regular_grid.resolution
    xmin, xmax, ymin, ymax, zmin, zmax = geo_model.solutions.grid.regular_grid.extent
    
    delx = (xmax - xmin)/nx
    dely = (ymax - ymin)/ny
    delz = (zmax - zmin)/nz
    
    # get unit IDs and restructure them
    #ids = np.round(geo_model.solutions.lith_block)
    #ids = ids.astype(int)
    ids = np.round(lithology_block)
    ids = ids.astype(int)
    
    liths = ids.reshape((nx, ny, nz))
    liths = liths.flatten('F')

    # group litho in space-saving way
    sequence = [len(list(group)) for key, group in it.groupby(liths)]
    unit_id = [key for key, group in it.groupby(liths)]
    combined = ["%s*%s" % (pair) for pair in zip(sequence,unit_id)]

    combined_string = " ".join(combined)
    
    # bcs
    with open(bcs_file, 'r') as file:
        bc_vals = file.read()
        file.seek(0)
        lines = len(file.readlines())
    
    # borehole logs
    if borehole_logs is not None:
        n_logs = len(borehole_logs)
        borehole_string = f"# borehole logs, records={n_logs} \n"
        for hole in range(n_logs):
            borehole_string += f"{borehole_logs[hole,0]}, {borehole_logs[hole,1]}, borehole{hole} \n"
    else:
        borehole_string = "!# borehole logs, records=0"
        
    if data_file is None:
        data_string = "!# data, records=0"
    else:
        with open(data_file, 'r') as file:
            data_vals = file.read()
            file.seek(0)
            data_lines = len(file.readlines())
        data_string = f"\n# data, records={data_lines-1}  !"
        data_string += data_vals
    # units
    unitstring = ""
    try:
        for index, rows in units.iterrows():
            unitstring += f"{rows['por']}    1.d0  1.d0  {rows['perm']}	 1.e-10  1.d0  1.d0  {rows['lz']}	0.  2077074.  10  2e-3	!{rows['surface']} \n" 	
    except:
        print("No units table found, filling in default values for petrophysical properties.")
        # get number of units and set units string
        units = geo_model.surfaces.df[['surface', 'id']]
        for index, rows in units.iterrows():
            unitstring += f"0.01d0    1.d0  1.d0  1.e-14	 1.e-10  1.d0  1.d0  3.74	0.  2077074.  10  2e-3	!{rows['surface']} \n" 	
        
    # input file as f-string
    fstring = f"""!==========>>>>> INFO
# Title
{filename}

# linfo
1 2 1 1

# runmode
1

# timestep control
0
1           1           0           0

# tunit
1
 
# time periods, records=1
0      60000000    200      lin
           
# output times, records=10
1
6000000
12000000
18000000
24000000
30000000
36000000
42000000
48000000
54000000
    
# file output: hdf

# active temp head

# PROPS=bas

# USER=none


# grid
{nx} {ny} {nz}

# delx
{nx}*{delx}

# dely
{ny}*{dely}

# delz
{nz}*{delz}

{borehole_string}

!==========>>>>> NONLINEAR SOLVER
# nlsolve
50 0

!==========>>>>> FLOW
# lsolvef (linear solver control)
1.d-8 64 500
# nliterf (nonlinear iteration control)
1.0d-6 1.0

!==========>>>>> TEMPERATURE
# lsolvet (linear solver control)
1.d-4 64 500
# nlitert (nonlinear iteration control)
1.0d-2 1.0

!==========>>>>> INITIAL VALUES
# temp init HDF5=temp_init.h5
! {nx*ny*nz}*11.0d0  

# head init HDF5=head_init.h5
! {nx*ny*nz}*7000

!==========>>>>> UNIT DESCRIPTION
!!
# units
{unitstring}

!==========>>>>>   define boundary properties
# temp bcd, simple=top, value=init

# temp bcn, simple=base, error=ignore
{nx*ny}*0.092

# head bcd, records={lines}
{bc_vals}

{data_string}

# uindex
{combined_string}"""

    if not path:
        path = './'
    if not os.path.exists(path):
        os.makedirs(path)

    f = open(path+filename, 'w+')
    
    f.write(fstring)
    f.close()
    
    print(f"Successfully exported geological model {filename} as SHEMAT-Suite input to "+path)      


# In[7]:


path_to_liths = '../../../../OneDrive/Dokumente/Work-GemPy/2021-05-04_GemPy/'


# In[8]:


lith_blocks = conc_lithblocks(path_to_liths)


# In[7]:


lith_blocks = np.load('../../../../OneDrive/Dokumente/Work-GemPy/2021-05-04_GemPy/concatenated/realistic_models_w_faults.npy')


# In[9]:


lith_blocks = np.load('../../../../OneDrive/Dokumente/Work-GemPy/2021-05-04_GemPy/base_model_w_fault.npy')


# In[10]:


np.unique(lith_blocks[0,:])


# In[8]:


lith_blocks.shape


# In[28]:


lith_blocks = np.load('../../../../data_var/2021-06-30_inference_miguel/inference_export/MetropolisAllStationsLithBlocks_topography.npy')


# In[11]:


lith_blocks_topo = np.array([])
for i in lith_blocks:
    lith_blocks_topo = np.append(lith_blocks_topo, topomask(geo_model, i))
lith_blocks_topo = lith_blocks_topo.reshape(len(lith_blocks), -1)


# In[12]:


lith_blocks_topo.shape


# In[14]:


np.save('../../../../data_var/2021-06-30_inference_miguel/inference_export/MetropolisAllStationsLithBlocks_topography.npy', lith_blocks_topo)


# In[13]:


np.save('../../../../OneDrive/Dokumente/Work-GemPy/2021-05-04_GemPy/concatenated/concatenated_w_topo.npy', lith_blocks_topo)


# In[21]:


xmin, xmax, ymin, ymax, zmin, zmax = geo_model.grid.regular_grid.extent


# In[22]:


head_bcd = '../../data/processed/GemPy/06_DTMs/Model_DTM_z_vals_SHEMAT_head_bc.txt'
temp_data = '../../data/processed/temperature_data/2021-05-11_temp_for_shemat_HRT_HT_BHT.csv'
temp_bcn = np.loadtxt('../../data/processed/temperature_data/basal_hf_gaussian')


# ## Set up the units for the SHEMAT-Suite model

# In[11]:


units = geo_model.surfaces.df[['surface', 'id']]
units


# In[12]:


tc = pd.read_csv('C:/Users/brigg/polybox/data_boreholes_aargau/interim/Loic/stat_TC_model_units_based_on_loic.csv', sep=';')
por = pd.read_csv('../../data/processed/petrophys/porosity_mean_median_stdev.csv')
perm = pd.read_csv('../../data/processed/petrophys/permeability_mean_median_stdev.csv')


# In[13]:


perm


# In[14]:


params = {'por': np.array([1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 0.0695, 0.0805, 0.0756, 0.05, 0.0198, 0.011, 0.115, 0.0625, 0.0074, 0.02, 0.01]).T,
         'perm': np.array([1e-16, 1e-16, 1e-16, 1e-16, 1e-16, 1e-16, 1.0e-16, 1.18e-17, 6.9e-14, 1.0e-15, 1.25e-15, 1.67e-16, 1.5e-15, 4.19e-14, 4.5e-15, 4.0e-15, 1e-18]),
         'lz':   np.array([2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.1, 1.94, 2.98, 1.93, 2.86, 1.95, 2.415, 2.98, 4.64, 2.87, 3.1])}


# In[15]:


units = units.join(pd.DataFrame(params, index=units.index))


# In[16]:


units


# In[17]:


air = {'surface': 'air',
       'id': units.shape[0]+1,
      'por': 1e-10,
      'perm': 1e-22,
      'lz': 100}
units = units.append(air, ignore_index=True)


# In[18]:


faults = {'surface': 'fault_units',
       'id': units.shape[0]+1,
      'por': 1e-1,
      'perm': 5e-14,
      'lz': 2.9}
units = units.append(faults, ignore_index=True)


# In[19]:


units


# In[24]:


shemade = ""
for c in range(len(lith_blocks)):
    model = lith_blocks[c,:]
    model_name = f"PCT_MC_{c}"
    export_shemat_suite_input_file(geo_model, lithology_block=model, units=units, 
                                   bcs_file=head_bcd, 
                                   data_file=temp_data,
                                   path='../../../../OneDrive/Dokumente/Work-GemPy/2021-05-04_GemPy/SHEMAT_MC_only_REALISTIC/',
                                  filename=model_name)
    shemade += model_name + " \n"
with open("../../../../OneDrive/Dokumente/Work-GemPy/2021-05-04_GemPy/SHEMAT_MC_only_REALISTIC/shemade.job", 'w') as jobfile:
    jobfile.write(shemade)


# In[21]:


export_shemat_suite_input_file(geo_model, lithology_block=lith_blocks, units=units, 
                               bcs_file=head_bcd,
                               data_file=temp_data,
                               path='../models/20210319_MC_no_middle_filling/',
                              filename='SHEMAT_PCT_base_model_fault')


# In[ ]:


model_shape = geo_model._grid.regular_grid.mask_topo.shape
topo_mask = geo_model._grid.regular_grid.mask_topo

ids = np.round(geo_model.solutions.lith_block)
ids = ids.astype(int)

liths = ids.reshape(geo_model._grid.regular_grid.mask_topo.shape)

liths[topo_mask] = np.unique(liths)[-1]+1

lith_blocks_topo = np.array([])

for i in lith_blocks:
    ids = np.round(i)
    ids = ids.astype(int)
    
    liths = ids.reshape(model_shape)
    liths[topo_mask] = np.unique(liths)[-1]+1
    lith_blocks_topo = np.append(lith_blocks_topo, liths)
lith_blocks_topo = lith_blocks_topo.reshape(n_iterations, -1)

