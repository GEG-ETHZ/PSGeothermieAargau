"""
Create a 3D model of a Permo-Carboniferous Trough (PCT)
=======================================================

Based on four seismic sections from the NAGRA report 
`NAGRA NTB 14-02 <https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NTBs\%202014\%20-\%202015/d_ntb14-02\%20Dossier\%20I.pdf>`_ [1], 
we extracted interface and orientation points of main eras (paleozoic, mesozoic, cenozoic) and major graben faults. 
Data from these 2D sections are complemented with data from GeoMol 2019, e.g. base of the PCT, thrusts, and normal faults. 

The lithological units comprise the permo-carboniferous filling (paleozoic), Mesozoic, Tertiary strata, as well as the crystalline basement rocks. An important decision before building the geological model,
is to define model units. Based on the purpose of the envisaged model, different units have to be defined. As the final result of this work will be an ensemble of advective heat-transport models,
key paremeters for defining units are permeability, porosity, thermal conductivity of different geological layers. As part of the exploration work of nagra 
(National Cooperative for the Disposal of Radioactive Waste), regional and local hydrogeological models were constructed. The therein defined hydrostratigraphy provides the basis for defining the 
model units of this geological model. The regional hydrogeologic model is presented in the report 
`NAGRA NAB 13-23 <https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NABs\%202004\%20-\%202015/e_nab13-023.pdf>`_ [2]. 

With the regional model covering an area comprising all potential storage sites defined by nagra, local models were built as well. These models comprise a more detailed hydrostratigraphy. 

The potential storage site "Jura Ost" is within our model area, thus we also consider the hydrostratigraphy defined in this local hydrogeological model presented in the report 
`NAGRA NAB 13-26 <https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NABs\%202004%20-\%202015/e_nab13-026.pdf>`_ [3].

The model comprises an area of 45 km x 32 km, in x- and y-direction, respectively. It extends down to a depth of 6 km, with reference sea level. 
This notebook demonstrates step-by-step how the model is generated within the open source modeling software `GemPy <https://www.gempy.org/>`_ [4].  
First, we will import libraries necessary to run this notebook:

"""

#%%

# Importing GemPy
import gempy as gp

# Import improved plotting features from GemPy
from gempy.plot import visualization_2d as vv
from gempy.plot import vista

# Importing auxilary libraries
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)

#%%
# This example code was generated with Gempy-Version:

print(f"GemPy Version: {gp.__version__}")

#%%
# Initialize the model
# --------------------
# For modeling the PermoCarboniferous trough (**PCT**) in GemPy, we need to initialize a GemPy model object. This model object comprises multiple input data, such as interface points and orientations, 
# which we previously stored in a `.csv` file. Further, we import the topography from a GeoTiff file.  
# Conceptually, we create two models:  
#
#.  1. With data of the the base of the PCT known  
#.  2. With additional data where the base of the PCT is inferred  
# 
# The distinction of inferred vs. known locations of the PCT is based on GeoMol 2019, an update geological model of the Swiss Molasse Basin and adjacent areas. Known and inferred parts of the PCT in 
# GeoMol can be seen `here <https://viewer.geomol.ch/webgui/gui2.php?viewHash=02171f57ee58a4082d3eb9cdc541c08b>`_.
# 
# In this notebook, the user can choose whether only the "known" parts of the PCT base will be considered for modeling, or also the the inferred parts.

#%%
# string either "known" or "inferred" to switch between model data
switch = "known" 

if switch == 'known':
    # Import data - NOT INFERRED
    # Create a model instance
    geo_model = gp.create_model('PCT_model')
    
    # Initialize the model, set dimension and load interface and orientation data
    gp.init_data(geo_model, [2640000, 2685000., 1240000., 1275000., -6000, 1000.], [50, 50, 50],
                path_i = '../../../Editorial-Transitional-Heatflow/data/processed/GemPy/00_gempy_inputs/2021-06-02_interfaces_no_fault_horizon_reduced_graben_and_mandach.csv',
                path_o = '../../../Editorial-Transitional-Heatflow/data/processed/GemPy/00_gempy_inputs/20201007_orientations_with_Jurathrust5_no_quat_meso_reduced2.csv')
    
    geo_model.set_topography(source='gdal', filepath='../../../Editorial-Transitional-Heatflow/data/processed/GemPy/06_DTMs/DTM_200_for_GemPy_Model.tif')
elif switch == 'inferred':
    # Import data - INFERRED
    # Create a model instance
    geo_model = gp.create_model('PCT_model_inferred')
    
    # Initialize the model, set dimension and load interface and orientation data
    gp.init_data(geo_model, [2640000, 2685000., 1240000., 1275000., -6000, 1000.], [50, 50, 50],
                path_i = '../../data/processed/GemPy/00_gempy_inputs/20201005_interfaces_Jurathrust5_pct_inferred.csv',
                path_o = '../../data/processed/GemPy/00_gempy_inputs/20201007_orientations_with_Jurathrust5_no_quat_meso_reduced2_pct_inferred.csv')
    
    geo_model.set_topography(source='gdal', filepath='../../data/processed/GemPy/06_DTMs/DTM_200_for_GemPy_Model.tif')


#%%
#
# To be coherent with existing geological models, e.g. geological cross-sections by nagra, we adapt the coloring for units according to 
# `NTB 14-02 <https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NTBs\%202014\%20-\%202015/d_ntb14-02\%20Dossier\%20I.pdf>`_ [5]. 
# For this, we create a color dictionary linking the units of the model to hex-color-codes.

col_dict = {'basement': '#efad83',
           'graben-fill': '#97ca68',
           'Mittlerer-Muschelkalk': '#f9ee3a',
           'Oberer-Muschelkalk': '#ffcf59',
           'Keuper': '#ffe19f',
           'Opalinuston': '#7f76b4',
           'Dogger': '#b0ac67',
           'Effinger-Schichten': '#47c4e2',
           'Malm': '#92d2ec',
           'USM': '#fbf379',
           'OMM': '#fbf379',
           'BIH-Basement-N': '#015482',
           'Fault-south': '#4585a8',
           'Fault_Basement_A': '#851515',
           'Vorwald_Basement': '#b54343',
           'Jurathrust5': '#5DA629',
           'Mandach': '#408f09'}

geo_model.surfaces.colors.change_colors(col_dict)

#%%
#
# Visualize the data distribution
# -------------------------------
# The following plot shows the different interface and orientation data loaded in the previous cell:

gp.plot_2d(geo_model, show_data=True, show_lith=False, show_results=False, direction='z', legend=False)

#%%
# The different colors in the plot represent the different model units. Circles represent the interface points, while arrows define the orientation of the respective surface in space. 
# 
# GemPy interpolates these input data in space using a universal co-kriging approach. Later on, we will set up the interpolator.

#%%
# Setting up Cross sections from the Nagra Report
# -----------------------------------------------
#
# As stated before, next to GeoMol [6], we incorporate geological interpretations from four migrated seismic sections, the NAGRA report 
# `NTB 14-02 <https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NTBs\%202014\%20-\%202015/d_ntb14-02\%20Dossier\%20I.pdf>`_. 
# For comparing the model results with the original interpretations, we define three cross sections in the model domain by specifying their start- and end-points and their resolution:
#
# set three sections which go roughly North South:

section_dict = {'section4_3':([2670826,1268793],[2676862,1255579],[100,100]),
                 'section4_4':([2649021,1267107],[2659842,1246715],[100,100]),
                 'section4_8':([2643284,1259358],[2680261,1268521],[100,100])}
geo_model.set_section_grid(section_dict)


#%%
# Display Model Information
# -------------------------
#
# In the following, we will go through model construction step-by-step. As an overview, we display the different units (here called `surfaces`) included in the model. 
# Note that also faults are surfaces within this model context. Currently, they are not marked as faults, and GemPy would treat them as the base of another geological model unit. 
# 
# To clarify, we model the base of a unit volume. That is, everything above the base surface is the respective unit, until the next base surface is reached. 
# In total, our model comprises 17 `surfaces`. Everything beneath is filled with the 18th surface, called `basement`. 
# 
# ### Surfaces
# The majority of the structural features, i.e. normal- and thrust faults, are named following the names in GeoMol.
# Main features of the model is the asymetric graben system, with the major normal faults (:code:`Fault_Basement_A`, :code:`Fault-south`, :code:`BIH-Basement-N`), 
# and the graben fill, which is not present beyond the graben shoulders, unless where it is inferred. 
# This, as well as the stop of major normal faults beneath the mesozoic units (the base of :code:`Mittlerer-Muschelkalk`) are important considerations for the modeling process. 

geo_model.surfaces


#%%
# Characteristics  
# ---------------
# One characteristic seen in the table above, is that all surfaces are assigned to a :code:`series` called :code:`Default series`. 
# A _series_ in GemPy indicates whether units should be interpolated using the same parameters. That is, all :code:`surfaces` within the same :code:`series` will be sub-parallel. 
# Thus, surfaces have to be grouped into different series, depending on their geometry in space. For instance, sub-parallel layers of a sedimentary sequence should be grouped in the same series, 
# while an unconformity, or a fault should be assorted to its own series. 
# 
# In this model, we group the majority of mesozoic and cenozoic units in one series, called :code:`Post_graben_series`. Only the mesozoic surface :code:`Mittlerer-Muschelkalk` will be assigned its own 
# series, as it forms the basal detachement of the Jura Mountains. Palaeozoic graben sediments are also assigned its own series.

# Assign formations to series
gp.map_series_to_surfaces(geo_model,
                         {"Thrust_Mandach": 'Mandach',
                          "Thrust_Jura": 'Jurathrust5',
                          #"Thrust_Jura6": 'Jurathrust6', #('Jurathrust4', 'Jurathrust5', 'Jurathrust6'),
                          "Fault_north_series": 'Fault_Basement_A',
                          "Fault_south_series": 'Fault-south',
                          "Vorwald_series": 'Vorwald_Basement',
                          "BIH_series": 'BIH-Basement-N',
                          "Fault_north_series": 'Fault_Basement_A',
                          "Fault_south_series": 'Fault-south',
                         "Post_graben_series": ('OMM',
                                                'USM',
                                                'Malm',
                                                'Effinger-Schichten',
                                                'Dogger', 
                                                'Opalinuston', 
                                                'Keuper',
                                                'Oberer-Muschelkalk'),
                          "Detachement": 'Mittlerer-Muschelkalk',
                         "Graben_series": 'graben-fill'},
                         remove_unused_series=True)
geo_model.surfaces

#%%
#
# Define Faults
# -------------
# To distinguish between lithological units and faults, we have to assign which series are faults. Faults can be infinite, i.e. have the same displacement throughout the model space, or they can be
# finite, meaning displacement will be less towards the fault edges (which are defined by the extent of interface points used as input).

geo_model.set_is_fault(['Thrust_Mandach', 'Thrust_Jura', 'Fault_north_series', 
                        'Fault_south_series', 'Vorwald_series', 'BIH_series'],
                      change_color=False)
geo_model.set_is_finite_fault(series_fault=['BIH_series', 'Vorwald_series'],
                              toggle=True)

#%%
# Bottom relation 
# ---------------
# To set whether a surface is eroding or not, we can set a series' `bottom_relation`. Per default, it is set to `Erosion`, meaning the base of a younger surface (higher up in the stratigraphic pile) 
# will cut through older surfaces. Setting the `bottom_relation` to `Onlap` will cause the opposite, i.e. younger surfaces stop on older ones.  
# We set the _Graben_series_ to onlap, as most of it is only present in the graben, i.e. hanging wall of the normal faults, but not in the foot wall.

geo_model.set_bottom_relation(series=['Graben_series'], bottom_relation='Onlap')

#%%
# Define Fault relations
# ----------------------
# With cross-cutting faults, we need to define fault relations, i.e. which fault stops at which. This is important, as some normal faults stop at others, e.g. :code:`BIH_Series` stops at 
# :code:`Fault_south_series`. Fault relations are set in a matrix, where :code:`True` sets that one fault stops at the other. If set to :code:`False` (the default), faults cross-cut each other 
# without any effects.
# 
# Further, fault relations are used to define whether a fault displaces lithological series, or not. Per default, all faults displace the lithological series, but not other faults. 
# This can be seen, if we plot the :code:`fault_relations` matrix:

geo_model.faults.faults_relations_df

#%%
# We know that faults do not affect all lithological series equally. For instance, thrusts will not affect the paleozoic sediments filling the graben. 
# Just as the mesozoic units are not affected by the normal faults. Thus we set up a fault relation matrix, considering:  
# ::
#    * thrusts only affect Mesozoic units    
#    * normal faults only affect Basement, Graben_series  
#    * normal faults stop at thrusts
# We can update the fault relations by creating a boolean matrix of shape similar to :code:`faults_relations_df`, to assign which fault displaces which unit, etc. Then we use this
# boolean matrix to set fault relations using the :code:`set_fault_relation()` method.

fr = np.array([[False, False, False, False, False, False, True,  False, False, False],
               [False, False, False, True,  False, False, True,  False, False, False],
               [False, False, False, False, True,  False,  False, True,  True, True],
               [False, False, False, False, False, True, False, False,  True, True],
               [False, False, False, False, False, False, False, True,  True, True],
               [False, False, False, False, False, False, False, False, True, True],
               [False, False, False, False, False, False, False, False, False, False],
               [False, False, False, False, False, False, False, False, False, False],
               [False, False, False, False, False, False, False, False, False, False],
               [False, False, False, False, False, False, False, False, False, False]])
geo_model.set_fault_relation(fr)

#%%
# Remember when we had a look at the input data and briefly mentioned the interpolator? We now set the interpolator function for the underlying co-kriging interpolation using theano:

gp.set_interpolator(geo_model,
                         compile_theano=True,
                         theano_optimizer='fast_compile',
                         verbose=[])

#%%
# Creating the model
# ------------------
# Now that we set the parameters and fault relations, it is time to start the modeling process. In Gempy, this is done using a single function :code:`gempy.comput_model` giving the prepared _geo_model_ 
# as input.

sol = gp.compute_model(geo_model, compute_mesh=True)

#%%
# For comparing model results with geological interpretations of the aforementioned seismic sections, we plot the model units on top of the seismic profiles. 
# Profiles 4.3 and 4.4 (nomenclature is taken from [1]) strike across the graben axis, while profile 4.8 goes roughly along the graben.
# 
# In the following plot, we model all profiles with the resulting geological grid, in the order from left to right: Profile 4.3, Profile 4.4, Profile 4.8.

gp.plot_2d(geo_model, section_names=list(section_dict), show_block=True, show_boundaries=False, show_data=False,
          show_topography=True, show_results=True)

#%%
# References
# ----------
# [1]: Naef, H., and Madritsch, H.: Tektonische Karte des Nordschweizer Permokarbontrogs: Aktualisierung basierend auf 2D-Seismik und Schweredaten. Nagra Arbeitsbericht NAB 14-017, (2014).  
# [2]: Gmünder, C., Malaguerra, F., Nusch, S., & Traber, D.: Regional Hydrogeo-logical Model of Northern Switzerland. Nagra Arbeitsbericht NAB, 13-23, (2014).  
# [3]: Luo, J., Monninkhoff, B., Becker J.K.: Hydrogeological model Jura Ost. Nagra Arbeitsbericht NAB, 13-26, (2014).  
# [4]: de la Varga, M., Schaaf, A., and Wellmann, F.: GemPy 1.0: Open-source stochastic geological modeling and inversion. Geoscientific Model Development, 12(1), (2019), 1. doi:http://dx.doi.org/10.5194/gmd-12-1-2019.  
# [5]: Gautschi, A., & Zuidema, P. (ed): Nagra technical report 14-02, geological basics-Dossier I-Introduction and summary; SGT Etappe 2: Vorschlag weiter zu untersuchender geologischer Standortgebiete mit zugehörigen Standortarealen für die Oberflächenanlage--Geologische Grundlagen--Dossier I--Einleitung und Zusammenfassung, (2014).  
# [6]: GeoMol Team (2015): GeoMol – Assessing subsurface potentials of the Alpine Foreland Basins for sustainable planning and use of natural resources – Project Report, 188 pp. (Augsburg, LfU).  