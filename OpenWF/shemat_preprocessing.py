#!/usr/bin/env python

"""
This file contains preprocessing methods for creating SHEMAT-Suite models, created using the 2-step-conditioning workflow developed in the project Pilot Study Geothermics Aargau.

Methods comprise masking the geological model with topography and calculation of parameters, such as the heatflow.
"""

# Libraries
import os,sys
import glob
import numpy as np
import itertools as it
import gempy as gp
import pandas as pd

__author__ = "Jan Niederau"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Jan Niederau"
__status__ = "Prototype"

def normalize(df):
    """Normalize a df-column

    Args:
        df ([type]): column of a pandas dataframe

    Returns:
        [type]: normalized version of said dataframe
    """
    return (df - df.min()) / (df.max() - df.min())


def topomask(geo_model, lith_block):
    """Method to mask the GemPy model with topography and change masked IDs to a new one for air.

    Args:
        geo_model (gp model): gempy geomodel
        lith_block (np.array): numpy array with the lith_block solution of a gempy model in a defined regular grid

    Returns:
        np.array: lith_block masked with topography, so including cell IDs representing topography. 
    """
    model_shape = geo_model._grid.regular_grid.mask_topo.shape
    topo_mask = geo_model._grid.regular_grid.mask_topo
    
    # round lithblock and change to integers
    ids = np.round(lith_block)
    ids = ids.astype(int)
    
    mask_lith = ids.reshape(model_shape)
    mask_lith[topo_mask] = np.unique(mask_lith)[-1] + 1
    
    return mask_lith

def conc_lithblocks(path: str='.'):
    """Concatenate multiple lith block files (npy files) in a folder into one numpy array

    Args:
        path (str, optional): Path to the numpy array files. Defaults to '.'.

    Returns:
        np.array: concatenated lith blocks in a numpy array
    """
    fids = glob.glob(path+'*.npy')
    lith_blocks = []
    
    for f in fids:
        blocks = np.load(f)
        lith_blocks.append(blocks)
        all_lith_blocks = np.concatenate(lith_blocks)
        
    return all_lith_blocks

def export_shemat_suite_input_file(geo_model, lithology_block, 
                                   units: pd.DataFrame=None, bcs_file: str=None, 
                                   data_file: str=None, borehole_logs: np.array=None,
                                   path: str=None, filename: str='geo_model_SHEMAT_input_erode'):
    """Method to export a 3D geological model as SHEMAT-Suite input-file for a conductive HT-simulation. 

    Args:
        geo_model (gp model): gempy model
        lithology_block (numpy array): array containing the lithology IDs for the regular grid of the model
        units (pd.DataFrame, optional): unit petrophysical parameters for SHEMAT-Suite model. Defaults to None.
        bcs_file (str, optional): boundary condition file for spatially varying boundary conditions (e.g. head by topography). Defaults to None.
        data_file (str, optional): data for calibrating the model, e.g. temperatures from boreholes. Defaults to None.
        borehole_logs (np.array, optional): coordinates for synthetic borehole logs, will write parameters such as pressure and temperature. Defaults to None.
        path (str, optional): save path for the SHEMAT-Suite input file. Defaults to None.
        filename (str, optional): name of the SHEMAt-Suite input file. Defaults to 'geo_model_SHEMAT_input_erode'.
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
    if bcs_file is not None:
        with open(bcs_file, 'r') as file:
            bc_vals = file.read()
            file.seek(0)
            lines = len(file.readlines())
            head_bcs = f"# head bcd, records={lines}\n{bc_vals}"
    else:
        head_bcs = f"# head bcd, simple=top, error=ignore\n{nx*ny}*{nz*delz}"
    
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
# temp init 
{nx*ny*nz}*11.0d0  

# head init 
{nx*ny*nz}*{nz*delz}

!==========>>>>> UNIT DESCRIPTION
!!
# units
{unitstring}

!==========>>>>>   define boundary properties
# temp bcd, simple=top, value=init

# temp bcn, simple=base, error=ignore
{nx*ny}*0.092

{head_bcs}

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