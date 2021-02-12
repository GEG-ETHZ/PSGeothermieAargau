# Methods for manipulating data files provided by the Canton Aargau
# developed by Katharina Kuhn und Jan Niederau
# 
# MIT License
# 
# Copyright (c) [2020] [Jan Niederau]
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture


def load_df(fid, header_lines=4, homogenize=True):

    dat = pd.read_csv(fid, header=header_lines, sep='\t')

    if fid[-1]=='C':
        new_col_names = ['Datum', 'Zeit', 'Temperatur(C)', 'Druck(mBar)', 'Tiefe(m)', 'down-uphole']
        dat.columns = new_col_names
    elif fid[-1]=='t':
        new_col_names = ['Tiefe(m)', 'Temperatur(C)']
        dat.columns = new_col_names

    return dat
    
def define_downhole_uphole(dat):
    #dat = dat.rename(columns={'Unnamed: 5': 'down-uphole'})
    
    dat['Tiefe(m)'] = dat['Tiefe(m)'].round(decimals=2)
    dat.loc[0,'down-uphole'] = 'down'

    for i in dat.index[1:]:
        if dat.loc[i,'Tiefe(m)'] >= dat.loc[i-1,'Tiefe(m)']:
            dat.loc[i,'down-uphole'] = 'down'
        else:
            dat.loc[i,'down-uphole'] = 'up'
    
    return dat

def reduce_df(dat, reduce_iters=3):
    # round depth to mm
    for i in range(reduce_iters):
        dat['Tiefe(m)'] = dat['Tiefe(m)'].round(decimals=2)
        dat_c = pd.DataFrame(columns=dat.columns)

        ## remove redundant points
        for i in dat.index[1:]:
            if dat['Tiefe(m)'].loc[i] > dat['Tiefe(m)'].loc[i-1]:
                dat_c = dat_c.append(dat.loc[i-1], ignore_index=True)
            else: #if dat['Tiefe(m)'].loc[i] <= dat['Tiefe(m)'].loc[i-1]:
                pass
    dat_c = dat_c.drop_duplicates(subset = ['Tiefe(m)'])

    return dat_c

def calc_gradT(dat):
    dz = np.gradient(dat['Tiefe(m)'][:])
    dt = np.gradient(dat['Temperatur(C)'][:])
    try:
        gradT = dt/dz
    except RuntimeWarning:
        print("divide by zero at {}".format(np.where(dz==0)))
    
    return gradT

def identify_peaks(dat, depthto=10, dist=1.):
    """Identify peaks in temperature over a certain depth interval

    Arguments:
        dat {dataframe} -- Dataframe with temperature and depth data

    Keyword Arguments:
        depthto {int} -- maximum depth of observation interval (in meters) (default: {10})
        dist {float} -- distance window of peak finding
    """
    temperature = dat.query("`Tiefe(m)` < @depthto")['Temperatur(C)']
    
    pks, _ = find_peaks(temperature, distance=dist)

    return pks

def cluster(dat, reduce_iters=3, ncluster=3, method='kmeans'):
    
    if reduce_iters != 0:
        dat = reduce_df(dat, reduce_iters)
    
    gradT = calc_gradT(dat)
    gradT[0] = gradT[1]
     
    x = np.stack((dat['Tiefe(m)'], gradT), axis=1)
    
    if method == 'kmeans':
        model = KMeans(n_clusters = ncluster)
        model.fit(x)
        y = model.predict(x)
    elif method == 'gmm':
        model = GaussianMixture(n_components = ncluster).fit(x)
        y = model.predict(x)

    dat['cluster'] = y
    dat['gradT'] = gradT
    
    return dat, model

def remove_zeros(dat):
    """Remove temperature gradients of 0

    Arguments:
        dat {dataframe} -- input dataframe, usually after calculating gradients

    Returns:
        [dataframe] -- cleaned dataframe
    """
    dat = dat[(dat != 0).all(1)].copy()
    
    gradT = calc_gradT(dat)
    dat['gradT'] = gradT

    return dat

def plot_temp(fid):
    dat = load_df(fid)
    label = str(fid).split('/')[-1]
    plt.plot(dat['Temperatur(C)'], -dat['Tiefe(m)'], label=label)
    plt.legend()
