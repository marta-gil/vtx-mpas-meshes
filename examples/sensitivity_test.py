import os

from vtxmpasmeshes.mesh_generator import full_generation_process
from vtxmpasmeshes.mpas_plots import compare_plot_mpas_regional_meshes, \
    view_mpas_regional_mesh


DATA_FOLDER = '/home/marta/PycharmProjects/vtx-mpas-meshes/data/' \
              'sensitivity-test-2'

# SITE: perdigao

lat_ref = 39.7136
lon_ref = -7.73

grids = []
for margin in [50, 75, 100, 125, 150, 250]:
    for size in [15, 25, 35, 50]:
        name = 'senst_s' + str(size).zfill(2) + '_m' + str(margin).zfill(3)
        folder = DATA_FOLDER + '/' + name + '/'
        basename = folder + name

        global_mesh = basename + '.grid.nc'
        regional_mesh = DATA_FOLDER + '/' + name + '.region.grid.nc'
        if not os.path.exists(regional_mesh):
            os.system('mkdir -p ' + folder)

            # I do a global mesh centered at 0,0 -> to use in MPAS workflow
            full_generation_process(
                global_mesh, 'doughnut',
                redo=False, do_plots=False, do_region=False,
                highresolution=3, lowresolution=20,
                size=size, margin=margin,
                lat_ref=0., lon_ref=0.,
            )

            # I do a regional mesh at my location -> with 4 layers
            full_generation_process(
                regional_mesh, 'doughnut',
                redo=False, do_plots=False, do_region=True,
                highresolution=3, lowresolution=20,
                num_boundary_layers=4,
                size=size, margin=margin,
                lat_ref=lat_ref, lon_ref=lon_ref,
            )

        f = DATA_FOLDER + '/' + name + '.resolution_mesh.png'
        if not os.path.isfile(f):
            print('Resolution Plots')
            view_mpas_regional_mesh(regional_mesh, vname='resolution',
                                    outfile=f)

        f = DATA_FOLDER + '/' + name + '.distortion_mesh.png'
        if not os.path.isfile(f):
            print('Distortion Plots')
            view_mpas_regional_mesh(regional_mesh, vname='cellDistortion',
                                    outfile=f, border_radius=size,
                                    cmap='magma')

        print('\n' + '*' * 30)
        print('\nDONE. This is the mesh ' + regional_mesh)
        grids.append(regional_mesh)


kwargs_set = {
    'default': {},
    '150': {'border_radius': 150, 'vmin': 0},
    'inner': {'border_radius': 30, 'vmin': 0},
}

for test, kwargs in kwargs_set.items():
    compare_plot_mpas_regional_meshes(grids,
                                      outfile=DATA_FOLDER +
                                              '/distortion.' + test + '.png',
                                      suptitle='Meshes comparison: Distortion',
                                      vname='cellDistortion',
                                      each_title='<NAME>: <NCELLS> cells',
                                      lat_ref=lat_ref,
                                      lon_ref=lon_ref,
                                      cmap='magma',
                                      **kwargs
                                      )

kwargs_set = {
    'default': {},
    '200': {'border_radius': 200, 'vmin': 3, 'vmax': 20},
    '50': {'border_radius': 50, 'vmin': 3, 'vmax': 15},
    'inner': {'border_radius': 30, 'vmin': 2.7, 'vmax': 3.3},
}

for test, kwargs in kwargs_set.items():
    compare_plot_mpas_regional_meshes(grids,
                                      outfile=DATA_FOLDER +
                                              '/compare.' + test + '.png',
                                      suptitle='Meshes comparison',
                                      each_title='<NAME>: <NCELLS> cells',
                                      lat_ref=lat_ref,
                                      lon_ref=lon_ref,
                                      **kwargs
                                      )
