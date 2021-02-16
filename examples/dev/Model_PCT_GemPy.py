"""
Export a GemPy Model
======================

Modeling the Permo-Carboniferous Trough in northern Switzerland, a potential geothermal reservoir system

"""

#%% 
# Create a 3D model of a Permo-Carboniferous Trough (PCT)
# =======================================================
#
# Based on four seismic sections from the NAGRA report [NAGRA NTB 14-02](https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NTBs\%202014\%20-\%202015/d_ntb14-02\%20Dossier\%20I.pdf) \cite{madritsch_nagra_2014}, we extracted interface and orientation points of main eras (paleozoic, mesozoic, cenozoic) and major graben faults. Data from these 2D sections are complemented with data from GeoMol 2019, e.g. base of the PCT, thrusts, and normal faults. 
# 
# The lithological units comprise the permo-carboniferous filling (paleozoic), Mesozoic, Tertiary strata, as well as the crystalline basement rocks. An important decision before building the geological model, is to define model units. Based on the purpose of the envisaged model, different units have to be defined. As the final result of this work will be an ensemble of advective heat-transport models, key paremeters for defining units are permeability, porosity, thermal conductivity of different geological layers. As part of the exploration work of nagra (National Cooperative for the Disposal of Radioactive Waste), regional and local hydrogeological models were constructed. The therein defined hydrostratigraphy provides the basis for defining the model units of this geological model. The regional hydrogeologic model is presented in the report 
# [NAGRA NAB 13-23](https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NABs\%202004\%20-\%202015/e_nab13-023.pdf) \cite{gmunder_regional_2014}. 
# 
# With the regional model covering an area comprising all potential storage sites defined by nagra, local models were built as well. These models comprise a more detailed hydrostratigraphy. 
# 
# The potential storage site "Jura Ost" is within our model area, thus we also consider the hydrostratigraphy defined in this local hydrogeological model presented in the report [NAGRA NAB 13-26](https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NABs\%202004%20-\%202015/e_nab13-026.pdf) \cite{luo_hydrogeological_2014}. 
# 
# The model comprises an area of 45 km x 32 km, in x- and y-direction, respectively. It extends down to a depth of 6 km, with reference sea level. This notebook demonstrates step-by-step how the model is generated within the open source modeling software [GemPy](https://www.gempy.org/) \cite{de_la_varga_gempy_2019}.  
# 
# First, we will import libraries necessary to run this notebook:

# In[1]:


# These two lines are necessary only if GemPy is not installed via pip
import sys, os
import numpy as np
sys.path = list(np.insert(sys.path, 0, r"../../../gempy/"))

# Importing GemPy
import gempy as gp

# Import improved plotting features from GemPy
from gempy.plot import visualization_2d as vv
from gempy.plot import vista

# Importing auxilary libraries
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib
matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)

# For embedding matplotlib figures
get_ipython().run_line_magic('matplotlib', 'inline')
import pyvista as pv
pv.set_plot_theme("document")


# # Initialize the model
# For modeling the PermoCarboniferous trough (**PCT**) in GemPy, we need to initialize a GemPy model object. This model object comprises multiple input data, such as interface points and orientations, which we previously stored in a `.csv` file. Further, we import the topography from a GeoTiff file.  
# Conceptually, we create two models:  
# 1. With data of the the base of the PCT known  
# 2. With additional data where the base of the PCT is inferred  
# 
# The distinction of inferred vs. known locations of the PCT is based on GeoMol 2019, an update geological model of the Swiss Molasse Basin and adjacent areas. Known and inferred parts of the PCT in GeoMol can be seen [here](https://viewer.geomol.ch/webgui/gui2.php?viewHash=02171f57ee58a4082d3eb9cdc541c08b).
# 
# In this notebook, the user can choose whether only the "known" parts of the PCT base will be considered for modeling, or also the the inferred parts.

# In[2]:


# string either "known" or "inferred" to switch between model data
switch = "known" 

if switch == 'known':
    # Import data - NOT INFERRED
    # Create a model instance
    geo_model = gp.create_model('PCT_model')
    
    # Initialize the model, set dimension and load interface and orientation data
    gp.init_data(geo_model, [2640000, 2685000., 1240000., 1275000., -6000, 1000.], [50, 50, 50],
                path_i = '../../data/processed/GemPy/00_gempy_inputs/20201005_interfaces_Jurathrust5_cleaned2_w_boreholes_fault_unit.csv',
                path_o = '../../data/processed/GemPy/00_gempy_inputs/20201007_orientations_with_Jurathrust5_no_quat_meso_reduced2.csv');
    #topo = geo_model.set_topography(source='gdal', filepath='../../data/processed/GemPy/Model_DTM_EPSG2056.tif');
    geo_model.set_topography(source='gdal', filepath='../../data/processed/GemPy/06_DTMs/DTM_200_for_GemPy_Model.tif');
elif switch == 'inferred':
    # Import data - INFERRED
    # Create a model instance
    geo_model = gp.create_model('PCT_model_inferred')
    
    # Initialize the model, set dimension and load interface and orientation data
    gp.init_data(geo_model, [2640000, 2685000., 1240000., 1275000., -6000, 1000.], [50, 50, 50],
                path_i = '../../data/processed/GemPy/00_gempy_inputs/20201005_interfaces_Jurathrust5_pct_inferred.csv',
                path_o = '../../data/processed/GemPy/00_gempy_inputs/20201007_orientations_with_Jurathrust5_no_quat_meso_reduced2_pct_inferred.csv');
    #topo = geo_model.set_topography(source='gdal', filepath='../../data/processed/GemPy/Model_DTM_EPSG2056.tif');
    geo_model.set_topography(source='gdal', filepath='../../data/processed/GemPy/06_DTMs/DTM_200_for_GemPy_Model.tif');


# In[3]:


col_dict = {'basement': '#efad83',
            'Graben-fill-low': '#07801a',
           'graben-fill': '#97ca68', ##bfdea1
           'Mittlerer-Muschelkalk': '#f9ee3a',
           'Oberer-Muschelkalk': '#ffcf59',
           'Keuper': '#ffe19f',
           'Opalinuston': '#7f76b4',
           'Dogger': '#b0ac67',
           'Effinger-Schichten': '#47c4e2',
           'Malm': '#92d2ec',
           'USM': '#fbf379',
           'OMM': '#fbf379',
           #'Basement_O_sued': '#017405',
           'BIH-Basement-N': '#015482',
           'Fault-south': '#4585a8',
           'Fault_Basement_A': '#851515',
           'Vorwald_Basement': '#b54343',
           'Jurathrust5': '#5DA629',
           'Mandach': '#408f09'}
geo_model.surfaces.colors.change_colors(cdict=col_dict)


# ## Visualize the data distribution
# The following plot shows the different interface and orientation data loaded in the previous cell:

# In[4]:


gp.plot_2d(geo_model, show_data=True, show_lith=False, show_results=False, direction='z', legend=False)


# The different colors in the plot represent the different model units. Circles represent the interface points, while arrows define the orientation of the respective surface in space. 
# 
# GemPy interpolates these input data in space using a universal co-kriging approach. As preparation, the interpolator has to be set up, which will be done in the following cell:

# In[5]:


# Set the interpolator function
# Create the theano model
gp.set_interpolator(geo_model,
                         compile_theano=True,
                         theano_optimizer='fast_compile',
                         verbose=[]);


# ## Setting up Cross sections from the Nagra Report
# 
# As stated before, next to GeoMol \cite{team_geomolassessing_2015}, we incorporate geological interpretations from four migrated seismic sections,  the NAGRA report [NAGRA NTB 14-02](https://www.nagra.ch/data/documents/database/dokumente/$default/Default\%20Folder/Publikationen/NTBs\%202014\%20-\%202015/d_ntb14-02\%20Dossier\%20I.pdf) \cite{madritsch_nagra_2014}. 
# For comparing the model results with the original interpretations, we define three cross sections in the model domain by specifying their start- and end-points and their resolution:

# In[6]:


# set 3 sections which go North South
section_dict = {'section4_3':([2670826,1268793],[2676862,1255579],[100,100]),
                 'section4_4':([2649021,1267107],[2659842,1246715],[100,100]),
                 'section4_8':([2643284,1259358],[2680261,1268521],[100,100])}
geo_model.set_section_grid(section_dict)


# ## Plot Data in 3D
# 
# Using PyVista, a high-level visualization software for 3D visualization, we can plot the data in 3D:

# In[7]:


gp.plot_3d(geo_model, show_surfaces=False, show_data=True, show_lith=False, show_topography=True,
           image=False, plotter_type='background', ve=1)


# # Display Model Information
# In the following, we will go through model construction step-by-step. As an overview, we display the different units (here called `surfaces`) included in the model. Note that also faults are surfaces within this model context. Currently, they are not marked as faults, and GemPy would treat them as the base of another geological model unit. 
# 
# To clarify, we model the base of a unit volume. That is, everything above the base surface is the respective unit, until the next base surface is reached. 
# In total, our model comprises 17 `surfaces`. Everything beneath is filled with the 18th surface, called `basement`. 
# 
# ### Surfaces
# The majority of the structural features, i.e. normal- and thrust faults, are named following the names in GeoMol.
# Main features of the model is the asymetric graben system, with the major normal faults (`Fault_Basement_A`, `Fault-south`, `BIH-Basement-N`), and the graben fill, which is not present beyond the graben shoulders, unless where it is inferred. This, as well as the stop of major normal faults beneath the mesozoic units (the base of `Mittlerer-Muschelkalk`) are important considerations for the modeling process. 

# In[8]:


geo_model.surfaces


# ## Characteristics  
# One characteristic seen in the table above, is that all surfaces are assigned to a `series` called `Default series`. A _series_ in GemPy indicates whether units should be interpolated using the same parameters. That is, all `surfaces` within the same `series` will be subparallel. Thus, surfaces have to be grouped into different series, depending on their geometry in space. For instance, sub-parallel layers of a sedimentary sequence should be grouped in the same series, while an unconformity, or a fault should be assorted to its own series. 
# 
# In this model, we group the majority of mesozoic and cenozoic units in one series, called `Post_graben_series`. Only the mesozoic surface `Mittlerer-Muschelkalk` will be assigned its own series, as it forms the basal detachement of the Jura Mountains. Palaeozoic graben sediments are also assigned its own series.

# In[9]:


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
                         "Graben_series": ('graben-fill', 'Graben-fill-low')},
                         remove_unused_series=True)
geo_model.surfaces


# ## Define Faults
# To distinguish between lithological units and faults, we have to assign which series are faults. 

# In[10]:


# Set Faults
geo_model.set_is_fault(['Thrust_Mandach', 'Thrust_Jura', 'Fault_north_series', 
                        'Fault_south_series', 'Vorwald_series', 'BIH_series'],
                      change_color=False);
geo_model.set_is_finite_fault(series_fault=['BIH_series', 'Vorwald_series'],
                              toggle=True);


# ### Bottom relation 
# To set whether a surface is eroding or not, we can set a series' `bottom_relation`. Per default, it is set to `Erosion`, meaning the base of a younger surface (higher up in the stratigraphic pile) will cut through older surfaces. Setting the `bottom_relation` to `Onlap` will cause the opposite, i.e. younger surfaces stop on older ones.  
# We set the _Graben_series_ to onlap, as most of it is only present in the graben, i.e. hanging wall of the normal faults, but not in the foot wall.

# In[11]:


geo_model.set_bottom_relation(series=['Graben_series'], bottom_relation='Onlap')


# ## Define Fault relations
# With cross-cutting faults, we need to define fault relations, i.e. which fault stops at which. This is important, as some normal faults stop at others, e.g. `BIH_Series` stops at `Fault_south_series`. Fault relations are set in a matrix, where `True` sets that one fault stops at the other. If set to `False` (the default), faults cross-cut each other without any effects.
# 
# Further, fault relations are used to define whether a fault displaces lithological series, or not. Per default, all faults displace the lithological series, but not other faults. This can be seen, if we plot the `fault_relations` matrix:

# In[12]:


# table of fault relations
geo_model.faults.faults_relations_df


# We know that faults do not affect all lithological series equally. For instance, thrusts will not affect the paleozoic sediments filling the graben. Just as the mesozoic units are not affected by the normal faults. Thus we set up a fault relation matrix, considering:  
# * thrusts only affect Mesozoic units    
# * normal faults only affect Basement, Graben_series  
# * normal faults stop at thrusts

# In[13]:


# with thrusts and mesozoic
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


# ## Creating the model
# Now that we set the parameters and fault relations, it is time to start the modeling process. In Gempy, this is done using a single function `gempy.comput_model` giving the prepared _geo_model_ as input.

# In[14]:


# Compute the model
sol = gp.compute_model(geo_model, compute_mesh=True)


# For comparing model results with geological interpretations of the aforementioned seismic sections, we plot the model units on top of the seismic profiles. Profiles 4.3 and 4.4 (nomenclature is taken from \cite{madritsch_nagra_2014}) strike across the graben axis, while profile 4.8 goes roughly along the graben.
# 
# ### Plot Profile 4.3

# In[15]:


p43 = vv.Plot2D(geo_model)
p43.create_figure((13,6));
sec_name = "section4_3"
a = p43.add_section(sec_name)

# Reading image
img = mpimg.imread('../../figs/Profil4-3_0-14einhalbkm.png')
# Plotting it inplace
a.imshow(img, origin='upper', alpha=.8, extent = (0, 14500, -2100, 1000))

p43.plot_topography(a, sec_name)
p43.plot_contacts(a, sec_name)


# ## Plot Profile 4.4

# In[16]:


p44 = vv.Plot2D(geo_model)
p44.create_figure((13,6));
sec_name = "section4_4"
a = p44.add_section(sec_name)

# Reading image
img = mpimg.imread('../../figs/Profil4-4_0-23km.png')
# Plotting it inplace
a.imshow(img, origin='upper', alpha=.8, extent = (0, 23000, -2100, 1000))

p44.plot_topography(a, sec_name)
p44.plot_contacts(a, sec_name)


# ## Plot Profile 4.8

# In[17]:


p48 = vv.Plot2D(geo_model)
p48.create_figure((13,6));
sec_name = "section4_8"
a = p48.add_section(sec_name)

# Reading image
img = mpimg.imread('../../figs/Profil4-8_0_38einhalbkm.png')
# Plotting it inplace
a.imshow(img, origin='upper', alpha=.8, extent = (0, 38500, -2100, 1000))

p48.plot_topography(a, sec_name)
p48.plot_contacts(a, sec_name)


# In comparison to the interpreted profiles, surface contours for faults and the graben sediments look promising. Mesozoic strata show folding instead of faulting by the thrusts. This has to be solved in the next iteration of the model. In Section 4.4, we see that `Jurathrust6` extends down into the basement (orange contour). While the surface extends down, it should be noted that `Jurathrust6` does not affect the paleozoic units, thus, having only a visualization effect. 
# In the next plot, we model all profiles with the resulting geological grid:

# In[18]:


gp.plot_2d(geo_model, section_names=list(section_dict), show_block=True, show_boundaries=False, show_data=False,
          show_topography=True, show_results=True)


# ## Plot 3D
# Of course, simulation results can also visualized in 3D. The following code-cells visualize the final block-model with topography, and as an interactive plot, where the user can slice through model planes:

# In[25]:


gp.plot_3d(geo_model, plotter_type='background', show_topography=True, show_surfaces=False, show_data=False,
          kwargs_plot_structured_grid={'opacity':1.0})


# In[26]:


gp.plot.plot_interactive_3d(geo_model, series=['Graben_series', 'Detachement'], show_topography=True)


# ## Gravity
# In addition to model geological structures, GemPy enables forward modeling of the gravity signal of the model. Forward gravity is calculated from the defined block model (here with resolution of 50 cells in each dimension) by applying the calculation method developed by Dezsö Nagy in 1966 \cite{nagy1966gravitational} for rectangular prisms in the vertical direction:  
# $$ F_z = G_{\rho}|||x\, ln(y + r) + y\, ln(x + r) - z\, arctan\bigg(\frac{xy}{zr}\bigg) |_{x_1}^{x_2}|_{y_1}^{y_2}|_{z_1}^{z_2}$$
# 
# Aside from density, the parameters stated in equation (1) are immutable. Hence, the decomposition of the according distance-related part of equation (1) ($t_z$) can be precalculated, yielding:
# 
# $$ F_z = G_{\rho} \cdot t_z $$.

# In[20]:


grav_res = 20
#[2635000, 2680000., 1240000., 1272000., -6000, 1000.]
X = np.linspace(2635000, 2680000., grav_res)
Y = np.linspace(1240000., 1272000., grav_res)
Z = 200.

xyz = np.meshgrid(X, Y, Z)
xy_ravel = np.vstack(list(map(np.ravel, xyz))).T


# For the forward simulation of gravity in this model, we create a custom, centered grid, where the majority of voxels are centered around a value, and resolution decreases towards the borders of the model / grid.  

# In[21]:


geo_model.set_centered_grid(xy_ravel, resolution = [10, 10, 15], radius=1000)


# Having set up a grid, the we need to assign densities to the different model units. Densities are taken from Abdelfettah et al. (2014) \cite{abdelfettah_characterization_2014} who did gravity simulation of Northern Switzerland, also including the PCT.

# In[22]:


# add densities - abdelfettah 2014
densities = [0, 0, 0, 0, 0, 0, 0, 2.5, 2.5, 2.5, 
             2.55, 2.55, 2.55, 2.55, 2.55, 2.6, 2.6, 2.57, 2.67]
geo_model.add_surface_values(densities, ['density'])


# In[23]:


gp.set_interpolator(geo_model, output=['gravity'], theano_optimizer='fast_run')


# In[24]:


sol = gp.compute_model(geo_model)
grav = sol.fw_gravity


# In[25]:


#gp._plot.plot_data(geo_model, direction='z')
plt.scatter(xy_ravel[:,0], xy_ravel[:,1], s=5)
im = plt.imshow(sol.fw_gravity.reshape(grav_res, grav_res), extent = (xy_ravel[:,0].min() + (xy_ravel[0, 0] - xy_ravel[1, 0])/2,
                                                       xy_ravel[:,0].max() - (xy_ravel[0, 0] - xy_ravel[1, 0])/2,
                                                       xy_ravel[:,1].min() + (xy_ravel[0, 1] - xy_ravel[30, 1])/2,
                                                       xy_ravel[:,1].max() - (xy_ravel[0, 1] - xy_ravel[30, 1])/2),
           cmap='viridis_r', origin='lower')
plt.colorbar(im, label='Bouguer gravity [mGal]')


# In[19]:


gp.plot_2d(geo_model, cell_number=42,
           direction='z', show_data=False, show_boundaries=False)


# ## Topology
# To distinguish between different conceptual models, we utilize the concept of model topology \cite{thiele_topology_2016-1}. Topological graphs of geological models are an abstract representation of model units which are in direct contact. Depending on the uniqueness of these graphs, we can categorize between different conceptual models.
# 
# For instance, the aforementioned input data of two conceptual models, with or without inferred data for the PCT base, would create models with different topological graphs. While model difference can also easily be estimated via visual inspection, the use of topological graphs facilitates analysis in case of many models. 
# 
# As this is the case in our study, where we analyse ensembles of uncertain geological structures, topological graphs allow for efficient ensemble analysis.

# In[20]:


from gempy.assets import topology as tp


# In[21]:


edges, centroids = tp.compute_topology(geo_model, cell_number=20, 
                                       direction='x', n_shift=1,
                                      voxel_threshold=5)


# In[23]:


gp.plot.plot_2d(geo_model, show_data=False, show_block=True, show_boundaries=False,
               show_topography=True, show_results=True, direction='x', cell_number=20)
gp.plot.plot_topology(geo_model, edges, centroids, direction='x')


# ## Save the Model

# In[21]:


geo_model.save_model(path='../../models/20200703_jordan_model_pct_inferred')


# # References
# 
# [<a id="cit-madritsch_nagra_2014" href="#call-madritsch_nagra_2014">1</a>] Madritsch H. and Deplazes G., ``_Nagra technical report 14-02, geological basics-Dossier II-Sediments and tectonic considerations; SGT Etappe 2: Vorschlag weiter zu untersuchender geologischer Standortgebiete mit zugehörigen Standortarealen für die Oberflächenanlage–Geologische Grundlagen–Dossier II–Sedimentologische und tektonische Verhältnisse_'', , vol. , number , pp. ,  2014.
# 
# [<a id="cit-gmunder_regional_2014" href="#call-gmunder_regional_2014">2</a>] Gmünder C., Malaguerra F., Nusch S. <em>et al.</em>, ``_Regional Hydrogeological Model of Northern Switzerland_'', , vol. , number 13-23, pp. ,  2014.  [online](https://www.nagra.ch/data/documents/database/dokumente/$default/Default%20Folder/Publikationen/NABs%202004%20-%202015/e_nab13-023.pdf)
# 
# [<a id="cit-luo_hydrogeological_2014" href="#call-luo_hydrogeological_2014">3</a>] Luo J., Monninkhoff B. and Becker J., ``_Hydrogeological model Jura Ost_'', , vol. , number 13-26, pp. ,  2014.  [online](https://www.nagra.ch/data/documents/database/dokumente/$default/Default%20Folder/Publikationen/NABs%202004%20-%202015/e_nab13-026.pdf)
# 
# [<a id="cit-de_la_varga_gempy_2019" href="#call-de_la_varga_gempy_2019">4</a>] de la Varga Miguel, Schaaf Alexander and Wellmann Florian, ``_GemPy 1.0: open-source stochastic geological modeling and inversion_'', Geoscientific Model Development, vol. , number , pp. ,  2019.
# 
# [<a id="cit-team_geomolassessing_2015" href="#call-team_geomolassessing_2015">5</a>] Team GeoMol, ``_GeoMol—assessing subsurface potentials of the Alpine Foreland Basins for sustainable planning and use of natural resources—project report_'', Augsburg: Bavarian Environment Agency, vol. , number , pp. ,  2015.
# 
# [<a id="cit-nagy1966gravitational" href="#call-nagy1966gravitational">6</a>] Nagy Dezs{\"o}, ``_The gravitational attraction of a right rectangular prism_'', Geophysics, vol. 31, number 2, pp. 362--371,  1966.
# 
# [<a id="cit-abdelfettah_characterization_2014" href="#call-abdelfettah_characterization_2014">7</a>] Abdelfettah Yassine, Schill Eva and Kuhn Pascal, ``_Characterization of geothermally relevant structures at the top of crystalline basement in Switzerland by filters and gravity forward modelling_'', Geophysical Journal International, vol. 199, number 1, pp. 226--241, oct 2014.  [online](https://academic.oup.com/gji/article/199/1/226/728762)
# 
# [<a id="cit-thiele_topology_2016-1" href="#call-thiele_topology_2016-1">8</a>] Thiele Samuel T., Jessell Mark W., Lindsay Mark <em>et al.</em>, ``_The topology of geology 1: Topological analysis_'', Journal of Structural Geology, vol. 91, number , pp. 27--38, October 2016.  [online](http://www.sciencedirect.com/science/article/pii/S0191814116301110)
# 
# 
