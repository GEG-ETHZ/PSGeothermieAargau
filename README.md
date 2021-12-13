# 📦 PSGeothermieAargau

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

> *Your documentation is a direct reflection of your software, so hold it to the same standards.*

This is the documentation repository for a sphinx gallery accompanying the [Pilot Study Geothermie Aargau](https://geg.ethz.ch/project-geothermal_aargau/).
The code written in this package comprises processing methods for an open-source workflow developed for this project. The project's aim is to study heat flow and subsurface temperature-distribution 
in a model-based stochastic approach. 
The gallery is meant to lead through the open workflow from data compilation to heat-transport simulation using step-by-step examples.  

## 🌟 Highlights

- Tutorial for leading through the developed workflow  
- Methods for creating input files for heat transport modeling in MOOSE or SHEMAT-Suite from GemPy Models
- Simple methods for calculating heat flow from temperature simulations
- Straight forward MonteCarlo rejection algorithm


## ℹ️ Overview

This repository 


### ✍️ Authors

The code and documentation was mainly developed by [Jan Niederau](https://github.com/Japhiolite) during his time at the [Geothermal Energy and Geofluids group](https://geg.ethz.ch)


## 🚀 Usage

```py
>>> import OpenWF.postprocessing as opp
>>> model_path = 'path/to/model/outputs/in/HDF5/format.h5'
>>> fid = pp.read_hdf_file(model_path, write=False)
>>> qx, qy, qz = opp.calc_cond_hf(fid)
print('calc_cond_hf yields conductive heat flow in all three model dimensions')
```


## ⬇️ Installation

The code project `OpenWF` is still in active development. Current procedure for installation is to clone this repository
and import it manually to your projects. A more convenient installation via pypi may come in the future.

#### Cloning the repository and adding it to the path

In a first step, clone this repository using the following command (or by manually downloading the zip file from the GitHub page)
(via HTTPS)
```bash
    $ git clone https://github.com/GEG-ETHZ/PSGeothermieAargau.git
```
(via SSH)
```bash
    $ git clone git@github.com:GEG-ETHZ/PSGeothermieAargau.git
```

using the python package `sys` you then append the path to the repository

```python    
    import sys
    sys.path.append("<path/to/cloned/repository/PSGeothermieAargau>")
```

followed by the import of the module

```python
    import OpenWF as owf
```

This package has been tested with Python 3.8 and works best currently with this pyhton Version.

`OpenWF` depends on a couple of other Python libraries:

* `GemPy` for creating the geological models
* `PyMC3` for bayesian inference (_optional_) Without PyMC3, only Massive Monte Carlo is possible.
* `matplotlib` for plotting
* `numpy` for efficient numerical methods
* `sqlite` for database access and manipulation  
* ...
  
an `environment.yml` file is provided for installation on Windows and Linux. These files can be used for creating
functioning anaconda environments for working with `OpenWF`.


## 💭 Feedback and Contributing

With regards to the scientific content, feel free to open a [disussion](https://github.com/GEG-ETHZ/PSGeothermieAargau/discussions). You are actively invited to open an [issue](https://github.com/GEG-ETHZ/PSGeothermieAargau/issues) for bugs/feature requests.

This project was meant as a pilot study to kick things off. Contributions to keep this project going and extend its capabilities is welcome. More detailed DEVELOPMENT and CONTRIBUTING guides will follow.


## 📜 License  

The content of this project documentation licensed under [CC-BY](https://choosealicense.com/licenses/cc-by-4.0/), the written source code of the developed workflow is licensed under [MIT](https://choosealicense.com/licenses/mit/).
