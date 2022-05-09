# vtx-mpas-meshes
Creation and analysis of MPAS-WRF meshes by Vortex

## Environment

    $ conda env create -n <envname> -f <path-to-environment.yml-file>

It contains:
* mpas-tools & jigsaw for mesh creation
* xarray, geopy, cartopy for mesh visualization
* packages for documentation
* jupyter notebook

The environment should be activated to run the scripts in this repository.

    $ conda activate <envname>

Then you can install the scripts of this repository by installing the `vtxmpasmeshes` source files using the `setup.py` file.

    (<envname>) $ pyhton setup.py install


## Mesh generation

Based on https://github.com/pedrospeixoto/MPAS-PXT/blob/master/grids/utilities/jigsaw/spherical_grid.py by Pedro S. Peixoto  ppeixoto@usp.br.

It uses:
* MPAS-Tools http://mpas-dev.github.io/MPAS-Tools/stable/mesh_creation.html#spherical-meshes
* jigsaw: https://github.com/dengwirda/jigsaw-python/tree/master/tests







