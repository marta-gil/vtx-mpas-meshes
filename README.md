# vtx-mpas-meshes
Creation and analysis of MPAS-WRF meshes by Vortex

## Environment

    $ conda env create -n <envname> -f <this-file.yml>

It contains:
* mpas-tools & jigsaw for mesh creation
* xarray, geopy, cartopy for mesh visualization
* packages for documentation
* jupyter notebook

The environment should be activated to run the scripts in this repository.

    $ conda activate <envname>

Then you can run, from the main directory of the repository, the scripts on any folder and have the results be saved on a folder called `data/`.

As an example (caution, the arguments may be outdated in newer versions):

    (<envname>) PATH/vtx-mpas-meshes$ python mesh-generation/mesh_generator.py -g doughnut -highr 3 -lowr 25 -size 30 -margin 100 -clat 42 -p -o
    (<envname>) PATH/vtx-mpas-meshes$ ls *
    environment.yml  README.md
    
    data:
    doughnut
    
    examples:
    
    mesh-generation:
    jigsaw_generator.py  mesh_generator.py  personalized_variable_resolution.py  __pycache__
    
    mesh-visualization:
    view_mpas_grid.py
    (<envname>) PATH/vtx-mpas-meshes$ ls data/doughnut/
    doughnut.graph.info     doughnut.grid.nc  
    doughnut-HFUN.msh       doughnut.jig  
    doughnut.log            doughnut-MESH.msh  
    doughnut.msh            doughnut.resolution.pdf  
    doughnut.triangles.nc


## Mesh generation

Based on https://github.com/pedrospeixoto/MPAS-PXT/blob/master/grids/utilities/jigsaw/spherical_grid.py by Pedro S. Peixoto  ppeixoto@usp.br.

It uses:
* MPAS-Tools http://mpas-dev.github.io/MPAS-Tools/stable/mesh_creation.html#spherical-meshes
* jigsaw: https://github.com/dengwirda/jigsaw-python/tree/master/tests







