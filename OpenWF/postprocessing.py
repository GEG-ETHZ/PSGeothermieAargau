#!/usr/bin/env python

"""
This file contains postprocessing methods for the SHEMAT-Suite models, created using the 2-step-conditioning workflow developed in the project Pilot Study Geothermics Aargau.

Methods comprise plotting and calculation of parameters, such as the heatflow.
"""

# import some libraries
from scipy.stats import hmean
import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_style('ticks')
sns.set_context('talk')


__author__ = "Jan Niederau"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Jan Niederau"
__status__ = "Prototype"

def read_hdf_file(filepath: str="."):
    """Short method to load an hdf5 file.

    Args:
        filepath (string): path to the stored .h5 file
    """
    model_fid = h5py.File(filepath, 'r')

    return model_fid

def available_parameters(file: h5py.File):
    """Return a summary of available parameters in the simulation file

    Args:
        file (h5py.File): loaded HDF5 file

    Returns:
        available parameters [dict]: dictionary with available parameters and short definition
    """

    all_params = {'comp': "compressibility",
    'delx': "discretization in x direction in meter",
    'dely': "discretization in y direction in meter",
    'delz': "discretization in z direction in meter",
    'df': "?",
    'ec': "?",
    'head': "hydraulic potential in meter",
    'itemp_bcd': "?",
    'itemp_bcn': "?",
    'kx': "log-permeability (square meter) in x direction",
    'ky': "log-permeability (square meter) in y direction",
    'kz': "log-permeability (square meter) in z direction",
    'lc': "?",
    'lx': "thermal conductivity in x direction in watt per meter and kelvin",
    'ly': "thermal conductivity in y direction in watt per meter and kelvin",
    'lz': "thermal conductivity in z direction in watt per meter and kelvin",
    'por': "porosity",
    'pres': "pressure in megapascal",
    'q': "?",
    'rc': "?",
    'rhof': "density water in kilogram per cubic meter",
    'temp': "temperature in degrees celsius",
    'temp_bcd': "temperature dirichlet boundary condition in degrees celsius",
    'temp_bcn': "temperature neumann boundary condition in degrees celsius",
    'uindex': "rock unit index - geological unit present in the cell",
    'visf': "fluid viscosity",
    'vx': "velocity in x direction in meters per second",
    'vy': "velocity in y direction in meters per second",
    'vz': "velocity in z direction in meters per second",
    'x': 'x coordinate in meters',
    'y': 'y coordinate in meters',
    'z': 'z coordinate in meters'}

    a_params = list(file.keys())

    a_subset = {key: value for key, value in all_params.items() if key in a_params}

    return a_subset

def plot_slice(file, parameter: str='temp', direction: str='x', cell_number: int=0, z_extent: float=6000.):
    """Plot a slice of available parameters through the model in x, y, or z direction. 

    Args:
        file (HDF5): hdf5 file of SHEMAT-Suite results
        parameter (str, optional): parameter to be plotted. . Defaults to 'temp'.
        direction (str, optional): [description]. Defaults to 'x'.
        cell_number (int, optional): [description]. Defaults to 0.
        z_extent (float, optional): [description]. Defaults to 6000..
    """
    if type(file)=='str':
        f = read_hdf_file(file)
    
    try:
        param = f[parameter][:,:,:]
    except KeyError:
        raise(f"Unable to open {parameter}. Does not exist in HDF5 File. See method 'available_parameters' for existing parameters.")

    uindex = f['uindex'][:,:,:]
    x = f['x'][0,0,:]
    y = f['y'][0,:,0]
    z = f['z'][:,0,0]
    zasl = z - z_extent
    
    if direction=='x':
        pa_cs = param[:,:,cell_number]
        ui_cs = uindex[:,:,cell_number]
        
        cs = plt.contourf(y,z-z_extent,pa_cs,23,cmap='viridis')
        plt.contour(y,z-z_extent,ui_cs, colors='#222222')
        plt.title(f'{parameter},x-direction, cell {cell_number}', fontsize=16)
        plt.tick_params(axis='both',labelsize=14)
        plt.xlabel('x [m]',fontsize=16)
        plt.ylabel('depth[m]',fontsize=16)
        cbar = plt.colorbar(cs,orientation='vertical')
        cbar.set_label(f'{parameter}',fontsize=16)
        cbar.ax.tick_params(labelsize=14)

        plt.show()
    
    elif direction=='y':
        pa_cs = param[:,cell_number,:]
        ui_cs = uindex[:,cell_number,:]
    
        cs = plt.contourf(x,z-z_extent,pa_cs,23,cmap='viridis')
        plt.contour(x,z-z_extent,ui_cs, colors='#222222')
        plt.title(f'{parameter},y-direction, cell {cell_number}', fontsize=16)
        plt.tick_params(axis='both',labelsize=14)
        plt.xlabel('y [m]',fontsize=16)
        plt.ylabel('depth [m]',fontsize=16)
        cbar = plt.colorbar(cs,orientation='vertical')
        cbar.set_label(f'{parameter}',fontsize=16)
        cbar.ax.tick_params(labelsize=14)

        plt.show()
        
    elif direction=='z':
        pa_cs = param[cell_number,:,:]
        ui_cs = uindex[cell_number,:,:]
        
        cs = plt.contourf(x,y,pa_cs,23,cmap='viridis')
        plt.contour(x,y,ui_cs, colors='#222222')
        plt.title(f'{parameter}, z-direction, {zasl[cell_number]} m depth', fontsize=16)
        plt.tick_params(axis='both',labelsize=14)
        plt.xlabel('x [m]',fontsize=16)
        plt.ylabel('y [m]',fontsize=16)
        cbar = plt.colorbar(cs,orientation='vertical')
        cbar.set_label(f'{parameter}',fontsize=16)
        cbar.ax.tick_params(labelsize=14)  

        plt.show() 

def find_nearest(array, value):
    """Find nearest cell index to given value

    Args:
        array ([type]): array with dimension values, i.e. x, y, or z direction 
        value ([type]): target value, i.e. some coordinate in x, y, or z direction

    Returns:
        [type]: closest cell index to defined value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    
    return idx

def calc_hf(data: h5py.File, depth_interval:list, model_depth: float=6000., direction: bool=False):
    """calculate the vertical heatflow over a certain depth interval.

    Args:
        data (h5py.File): [description]
        depth_interval (list): list of depth interval with [deeper, shallower] values
        model_depth (float, optional): vertical extent of the model in meter. Defaults to 6000..
        direction (boolean, optional): if set to True, direction of heatflow will be included, i.e. negative heat flows for outward ones

    Returns:
        hf [array]: average heat flow over the defined depth interval
    """
    zasl = data['z'][:,0,0] - model_depth
    upper = find_nearest(zasl, depth_interval[0])
    lower = find_nearest(zasl, depth_interval[1])
    
    #temp_diff = data['temp'][upper,:,:] - data['temp'][lower,:,:]
    temp_diff = np.sum(np.gradient(data['temp'][upper:lower+1,:,:], axis=0),axis=0)
    tc_av = hmean(data['lz'][upper:lower+1,:,:])
    #z_diff = zasl[upper] - zasl[lower]
    z_diff = np.sum(np.gradient(zasl[upper:lower+1], axis=0),axis=0)
    if direction==True:
        hf = - tc_av * (temp_diff/z_diff)
    else:
        hf = np.abs(-tc_av * (temp_diff/z_diff))

    return hf

def plot_logs(data, delz, borehole: int=0, apa=0.7, z_extent: float=6000., col_dict: dict={}):
    """Plot temperature logs of simulated and observed temperatures with lithologies as background

    Args:
        data (pandas Dataframe)): dataframe with borehole temperature data and cell values, i.e. i,j,k
        delz (float): cell size in z direction
        borehole (int, optional): borehole ID. Defaults to 0.
        apa (float, optional): alpha value. Defaults to 0.7.
        z_extent (float, optional): vertical model dimension. Defaults to 6000..
    """
    col_dict = {20: '#efad00',
            18: '#efad83',
            17: '#609133',
           16: '#97ca68', 
           15: '#f9ee3a',
           14: '#ffcf59',
           13: '#ffe19f',
           12: '#7f76b4',
           11: '#b0ac67',
           10: '#47c4e2',
           9: '#92d2ec',
           8: '#fbf379',
           7: '#fbf379',
           19: '#fffafa'}
    
    unique_boreholes = data[['i','j']].drop_duplicates().reset_index(inplace=False)
    indx = unique_boreholes.iloc[borehole]
    
    get_hole_data = data.query(f"i=={indx['i']} and j=={indx['j']}")
    litho_changes = np.where(get_hole_data['unit'].values[:-1] != get_hole_data['unit'].values[1:])[0]
    litho_f_color = list(get_hole_data['unit'].values[litho_changes])
    litho_f_color.append(get_hole_data['unit'].iloc[-1])
    
    depth = (get_hole_data['k'].values*delz) - (z_extent - delz/2)
    plt.scatter(get_hole_data['calc'], depth, s=80, marker='o', edgecolor='gray', facecolor='#055ff4')
    plt.scatter(get_hole_data['obs'], depth, s=60, marker='+', edgecolor='black', facecolor='#faa00b')
    
    for l in range(len(depth[litho_changes])):
        try:
            plt.axhspan(depth[litho_changes][l+1], 
                       depth[litho_changes][l], 
                        facecolor=col_dict[litho_f_color[l]], alpha=apa, zorder=-1)
        except IndexError:
            plt.axhspan(depth[0], depth[litho_changes][0],
                       facecolor=col_dict[litho_f_color[0]], alpha=apa, zorder=-1)
            plt.axhspan(depth[-1], depth[litho_changes][l],
                       facecolor=col_dict[litho_f_color[-1]], alpha=apa, zorder=-1)
            
def load_inv(filepath: str='.'):
    """Load SHEMAT-Suite simulation data files

    Args:
        filepath (str, optional): path to shemat-suite datafile. Files include simulated and observed values for #data nodes in the shemat-suite model. Defaults to '.'.

    Returns:
        [Pandas Dataframe]: Dataframe with loaded values
    """
    data = pd.read_csv(filepath, skipinitialspace=True, skiprows=2, header=[0], sep=" ")
    cols = list(data.columns.values[1:])

    data = data.dropna(axis=1, how='all')
    data.columns = cols
    
    return data

def fahrenheit_to_celsius(temp_fahrenheit, difference=False):
    """Transform Fahrenheit temperatures to Celsius

    Args:
        temp_fahrenheit (float): temperature in fahrenheit
        difference (bool, optional): relative difference to zero. Defaults to False.

    Returns:
        [type]: temperature in Celsius
    """
    if not difference:
        return (temp_fahrenheit - 32) * 5 / 9
    else:
        return temp_fahrenheit * 5 / 9
      
def read2dict(file: str='.'):
    """Read a parameter file and store it as a dictionary

    Args:
        file (string): path to file loaded. These are the parameter files generated by SHEMAT-Suite in gradient based inversion

    Returns:
        data_dict [dictionary]: dictionary with parameter blocks as arrays for different sections in the parameter file.
    """
    data_dict = {}
    with open(file) as f:
        for i, group in enumerate(get_groups(f, "#"), start=1):
            data_dict[group[0]]=group[1:]
            
    for name in list(data_dict.keys()):    
        for i in range(len(data_dict[name])):
            str_list = data_dict[name][i].split(' ')
            str_list = list(filter(None, str_list))
            unit_idx = str_list.index('unit')
        
            del str_list[unit_idx:]
            data_dict[name][i] = np.array(str_list, dtype=float)
        
        data_dict[name] = np.array(data_dict[name])
                
    return data_dict

def plot_log_inv(inv_data, forw_data, parameter_data, borehole, delz, z_extent: float, col_dict: dict, reduction: int=1, apa=0.7, save=False):
    """Plot the simulated and measured temperature log of a borehole in front of borehole stratigraphy

    Args:
        inv_data (pandas dataframe): dataframe with simulated data from shemat-suite
        forw_data (pandas dataframe): dataframe with observed data
        parameter_data (dictionary): dictionary with arrays from parameter files, read in with read2dict method
        borehole (int): borehole ID
        delz (float): vertical cell size of the model in meters (assuming regular grid)
        z_extent(float): vertical model extent in meters
        col_dict (dictionary): dictionary pairing lithology IDs to colors
        reduction (int): data reduction for plotting (in case you have so many points in your log)
        apa (float, optional): alpha value. Defaults to 0.7.
        save (bool, optional): save option for plot. Defaults to False.
    """
      
    bhole_data = forw_data.query("Borehole_Name==@borehole")
    depth = (inv_data['k']*delz) - z_extent
    litho_changes = np.where(inv_data['unit'].values[:-1] != inv_data['unit'].values[1:])[0]
    
    litho_f_color = list(inv_data['unit'][litho_changes].values)
    litho_f_color.append(inv_data['unit'].values[-1])
    
    fig = plt.figure(figsize=[10,8])
    apa = 0.7
    plt.scatter(bhole_data['Temperatur[C]'][::reduction], bhole_data['Z[asl]'][::reduction], 
                marker='^', facecolor='white', edgecolors='black', zorder=0)
    
    plt.plot(inv_data['calc'], depth, '--', c='#FF1111', linewidth=3, zorder=1)
    
    #plt.hlines(depth[litho_changes],0,100, color='black')
    for i in range(len(depth[litho_changes].values)):
        try:
            plt.axhspan(depth[litho_changes].values[i+1], 
                       depth[litho_changes].values[i], 
                        facecolor=col_dict[litho_f_color[i+1]], alpha=apa, zorder=-1)
        except IndexError:
            plt.axhspan(depth.iloc[0], depth[litho_changes].values[0],
                       facecolor=col_dict[litho_f_color[0]], alpha=apa, zorder=-1)
            plt.axhspan(depth.iloc[-1], depth[litho_changes].values[i],
                       facecolor=col_dict[litho_f_color[-1]], alpha=apa, zorder=-1)
            
    plt.text(bhole_data['Temperatur[C]'].min(),depth.values[-1]+50,'Calc. HF [mW/m$^2$] = {:.3f}'.format(parameter_data['# bcunits (aposteriori), iter =   final\n'][0][1]))
    
    plt.xlabel('Temperature [°C]')
    plt.ylabel('Depth [m asl.]')
    
    if save:
        fig.savefig(f"{borehole}_inv_plot_every_{reduction}_point.png", dpi=300, bbox_inches='tight')
        
def plot_apriori_aposteriori(parameter_data, parameter, borehole, save=False, interval=[0,100]):
    
    param = param_dict[parameter]
    
    y = [2, 1]
    my_ticks = ['apriori', 'aposteriori']
    
    apriori = parameter_data['# apriori (units apriori)\n'][:, param]
    aposteriori = parameter_data['# units (aposteriori), iter =   final\n'][:, param]
    
    for i in range(len(apriori)):
        try:
            plt.scatter(apriori[i], y[0], c=col_dict[i+1], edgecolors='black')
            plt.scatter(aposteriori[i], y[1], c=col_dict[i+1], edgecolors='black')
        except KeyError:
            plt.scatter(apriori[i], y[0], c=col_dict[i], edgecolors='black')
            plt.scatter(aposteriori[i], y[1], c=col_dict[i], edgecolors='black')
    
    plt.yticks(y, my_ticks)
    plt.xlabel(f"{parameter}")
    plt.xlim(interval)
    
    if save:
        fig.savefig(f"{borehole}_parameter_{parameter}.png", dpi=300, bbox_inches='tight')
  
def extract_parameters(datafile, parameters: list=['temp','uindex'], dimension: int=3, direction: str='x'):
    """extract single parameter fields from an HDF5 file for postprocessing.

    Args:
        datafile ([type]): [description]
        dimension (int, optional): [description]. Defaults to 3.
        direction (str, optional): [description]. Defaults to 'x'.

    Returns:
        [type]: [description]
    """
    param_dict = {}
    if type(datafile)==str:
        f = read_hdf_file(datafile)
    else:
        f = datafile
    
    z, y, x = f['temp'].shape
    param_dict['x'] = f['x'][0,0,:]
    param_dict['y'] = f['y'][0,:,0]
    param_dict['z'] = f['z'][:,0,0]
    
    
    if dimension==3:
        for i in parameters:
            param_dict[i] = f[i][:,:,:]
    elif dimension==2:
        if direction=='x':
            for i in parameters:
                param_dict[i] = f[i][:,:,x//2]
        elif direction=='y':
            for i in parameters:
                param_dict[i] = f[i][:,y//2,:]
        elif direction=='z':
            for i in parameters:
                param_dict[i] = f[i][z//2,:,:]
    return param_dict

# def homogenize_comment(file):
#     inp = open(file, 'rt')
#     outp = open(file+'_homogenized', 'wt')
#     for line in inp:    
#         outp.write(line.replace('%', '#'))#.replace('    ',','))
        
#     inp.close()
#     outp.close()            

# def get_groups(seq, group_by):
#     data = []
#     for line in seq:
#         # Here the `startswith()` logic can be replaced with other
#         # condition(s) depending on the requirement.
#         newline = line.replace('%', '#')
            
#         if newline.startswith(group_by):
#             if data:
#                 yield data
#                 data = []
#         data.append(newline)

#     if data:
#         yield data