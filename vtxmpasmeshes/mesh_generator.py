import os
import time

import xarray as xr
import numpy as np

from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.mesh.conversion import convert
from mpas_tools.io import write_netcdf

from vtxmpasmeshes.jigsaw_generator import jigsaw_gen_sph_grid
from vtxmpasmeshes.plot_utilities import view_resolution_map, \
    plot_mpas_darray
from vtxmpasmeshes.dataset_utilities import distance_latlon_matrix, \
    open_mpas_regional_file

PATH_LIMITED_AREA = '/home/marta/PycharmProjects/MPAS-Limited-Area'


def doughnut_variable_resolution(distances, highresolution,
                                 lowresolution, size, margin, **kwargs):

    # initialize with a very low resolution
    final_res_dist = 1000
    resolution = final_res_dist * np.ones(distances.shape)

    # Inner Circle
    # ------------------------------
    # inner circle of radius <size>km and resolution <highresolution>km
    # ie, where distances_array < size set resolution = highresolution
    central_zone = distances <= size
    resolution = np.where(central_zone, highresolution, resolution)

    # 1st Variable Resolution Ring
    # -------------------------------
    # Starts at the edge of the inner circle and lasts for <margin>km,
    # so it ends at a distance radius = <size>+<margin> km
    # The resolution in this area transitions linearly from
    # <highresolution> to <lowresolution>.

    # distance where the 1st variable resolution ring finishes
    radius = size + margin
    # slope at the 1st variable resolution ring
    slope = (lowresolution - highresolution) / margin

    transition_zone = (distances > size) & (distances <= radius)
    distance_to_innercircle = distances - size
    transition_vals = highresolution + slope * distance_to_innercircle
    resolution = np.where(transition_zone, transition_vals, resolution)

    # Fixed Low Resolution ring
    # -------------------------------
    # Starts at the edge of the 1st variable resolution ring (<radius>)
    # circle and lasts for "10 rings of cells", which we can approximate
    # to buffer=10*<lowresolution>km, so it ends at a distance
    # border = <radius>+<buffer> km
    # The resolution in this area is constant: <lowresolution> km

    # width of the low resolution ring
    buffer = 10*lowresolution
    # distance where the low resolution ring finishes
    border = radius + buffer

    lowres_zone = (distances > radius) & (distances <= border)
    resolution = np.where(lowres_zone, lowresolution, resolution)

    # Crazy increase in cell size ring
    # -------------------------------
    # Further than the <border> of the low resolution area we
    # keep increasing the resolution (using the same <slope>) as
    # before. This is just done to avoid having a very big global mesh.
    # The regional mesh generated after cutting this global mesh
    # will not include this part.
    # We stop this ring when the variable resolution reaches the
    # maximum width cell final_res_dist (which is the value we
    # initialized the resolution array with).

    distance_to_lowring= distances - border
    second_transition_vals = lowresolution + slope * distance_to_lowring

    second_transition_zone = (distances > border) & \
                             (second_transition_vals < final_res_dist)

    resolution = np.where(second_transition_zone,
                          second_transition_vals, resolution)

    # Finished!

    return resolution


def variable_resolution_latlonmap(grid, **kwargs):

    print('\n>> Creating a variable resolution map')

    # GRID
    # Create a global lat/lon grid at high resolution

    highresolution = kwargs.get('highresolution', 10.)  # grid size in km
    print('\tResolution in km of lat/lon grid: %.1f' % highresolution)

    dist_degrees = highresolution / 110.

    nlat = int(180. / dist_degrees) + 1
    nlon = int(360. / dist_degrees) + 1

    ds = xr.Dataset(
        coords={
            'lat': np.linspace(-90., 90., nlat),
            'lon': np.linspace(-180., 180., nlon),
        }
    )

    # DISTANCE
    # Compute distance from each point to the reference point

    lat_ref = kwargs.get('lat_ref', None)
    if lat_ref is None:
        lat_ref = 0.
    lon_ref = kwargs.get('lon_ref', None)
    if lon_ref is None:
        lon_ref = 0.
    kwargs.update({'lat_ref': lat_ref, 'lon_ref': lon_ref})

    print('\tComputing the distance to the reference point '
          '(%.2f, %.2f)' % (lat_ref, lon_ref))
    dists = distance_latlon_matrix(ds.coords['lat'], ds.coords['lon'],
                                   lat_ref=lat_ref, lon_ref=lon_ref,
                                   do_tile=True)

    ds['distance'] = xr.DataArray(data=dists, dims=('lat', 'lon'))

    # RESOLUTION
    # Set the resolution value at each point
    print('\tComputing resolutions using technique %s' % grid)

    if grid == 'doughnut':

        # Setting parameters to the defaults if not passed as arguments
        defaults = {'lowresolution': 3,
                    'highresolution': 25,
                    'size': 40,
                    'margin': 20}

        for name, default in defaults.items():
            kwag = kwargs.get(name, None)
            if kwag is None:
                kwargs.update({name: default})

        kwargs['radius'] = kwargs['size'] + kwargs['margin']
        kwargs['buffer'] = 10*kwargs['lowresolution']
        kwargs['border'] = kwargs['radius'] + kwargs['buffer']

        print('\t' + repr(kwargs))
        res = doughnut_variable_resolution(ds['distance'].values,
                                           **kwargs)
        ds['resolution'] = xr.DataArray(data=res, dims=('lat', 'lon'))

    else:
        raise ValueError('!! Grid %s not implemented.' % grid)

    ds.attrs = kwargs

    return ds


def get_mesh_from_resolution(resolution_ds, basename='./mesh'):
    print('\n>> Generating an MPAS mesh')

    # jigsaw
    print('\n\t .- Jigsaw generation')
    mesh_file = jigsaw_gen_sph_grid(resolution_ds['resolution'].values,
                                    resolution_ds['lon'].values,
                                    resolution_ds['lat'].values,
                                    basename=basename)

    # mpas-tools

    print('\n\t .- Jigsaw to netCDF')
    out_file_triangles = basename + '.triangles.nc'
    jigsaw_to_netcdf(msh_filename=mesh_file,
                     output_name=out_file_triangles,
                     on_sphere=True, sphere_radius=1.0)

    print('\n\t .- Convert to MPAS format')
    out_file_mpas = basename + '.grid.nc'
    write_netcdf(
        convert(xr.open_dataset(out_file_triangles),
                dir=os.path.dirname(basename),
                graphInfoFileName=basename + ".graph.info"),
        out_file_mpas)

    return out_file_mpas


def cut_circular_region(mpas_global_file,
                        region_radius_meters,
                        regional_grid=None,
                        regional_grid_info=None,
                        lat_cen=0., lon_cen=0.,
                        path_create_region=PATH_LIMITED_AREA,
                        ):

    print('\n>> Cutting a circular region')
    print('\t centered at %.4f, %.4f' % (lat_cen, lon_cen))
    print('\t with radius %.1fkm' % float(region_radius_meters/1000))

    if region_radius_meters < 5000:
        raise ValueError('Do you want a %.0fm radius region?'
                         'That looks way too small. Maybe you passed '
                         'the radius in km instead as in meters.'
                         % region_radius_meters)

    if not os.path.exists(mpas_global_file):
        raise IOError('Wanted to use the MPAS global file %s but'
                      'it does not seem to exist.' % mpas_global_file)

    os.system('cp ' + mpas_global_file + ' global.grid.nc')

    with open('points.txt', 'w') as f:
        f.write('Name: circle\n')
        f.write('Type: circle\n')
        f.write('Point: %.4f, %.4f\n' % (lat_cen, lon_cen))
        f.write('radius: %.0f\n' % region_radius_meters)

    # If there are regional files we should erase them
    os.system('rm -f circle.grid.nc')
    os.system('rm -f circle.graph.info')

    # Execute create region
    os.system(path_create_region + '/create_region points.txt ' +
              mpas_global_file)

    if not os.path.exists('circle.grid.nc') or \
            not os.path.exists('circle.graph.info'):
        raise AttributeError('The regions were not generated correctly')

    if regional_grid is not None:
        os.system('cp circle.grid.nc ' + regional_grid)

        if regional_grid_info is None:
            regional_grid_info = regional_grid.replace('.nc', '.graph.info')
    else:
        regional_grid = 'circle.grid.nc'

    if regional_grid_info is not None:
        os.system('cp circle.graph.info ' + regional_grid_info)
    else:
        regional_grid_info = 'circle.graph.info'

    return regional_grid, regional_grid_info


def full_generation_process(mpas_grid_file, grid, redo=True,
                            do_plots=True, **kwargs):

    if os.path.isfile(mpas_grid_file) and not redo:
        print(' .. already available')
        return

    graph_info_file = mpas_grid_file.replace('.nc', '.graph.info')
    path_save = os.path.dirname(mpas_grid_file)

    os.system('rm -f ' + mpas_grid_file)
    os.system('rm -f ' + graph_info_file)

    start_time = time.time()
    resolution_ds = variable_resolution_latlonmap(grid, **kwargs)
    duration_resolution = time.time() - start_time
    print(' .. finished finding resolution map: %.3fs\n\n' % duration_resolution)

    border = resolution_ds.attrs['border']
    radius = resolution_ds.attrs['radius']
    lat_ref = resolution_ds.attrs['lat_ref']
    lon_ref = resolution_ds.attrs['lon_ref']

    if do_plots:
        print('Plotting')
        start_time = time.time()
        view_resolution_map(resolution_ds,
                            pdfname=path_save + '/resolution.pdf',
                            list_distances=[
                                #1000,
                                500, border, radius])
        duration_plots = time.time() - start_time
        print(' .. finished doing resolution plots: %.3fs\n\n' % duration_plots)

    start_time = time.time()
    tmp_mesh_file = get_mesh_from_resolution(resolution_ds,
                                             basename='mesh')
    duration_gen = time.time() - start_time
    print(' .. finished generating global mesh: %.3fs\n\n' % duration_gen)

    mpas_grid_file_tmp = mpas_grid_file + '.tmp'
    start_time = time.time()
    cut_circular_region(tmp_mesh_file, radius*1000,
                        regional_grid=mpas_grid_file_tmp,
                        regional_grid_info=graph_info_file,
                        lat_cen=lat_ref, lon_cen=lon_ref)
    duration_region = time.time() - start_time
    print(' .. finished cutting region: %.3fs\n\n' % duration_region)

    if not os.path.isfile(mpas_grid_file_tmp) or \
            not os.path.isfile(graph_info_file):
        raise IOError('The file we had to generate was not generated')

    # Open dataset and update attributes
    mpas_ds = xr.open_dataset(mpas_grid_file_tmp)
    for name, value in resolution_ds.attrs.items():
        mpas_ds.attrs['vtx-param-' + str(name)] = value

    durations_process = {
        'resolution': duration_resolution,
        'generation': duration_gen,
        'region': duration_region,
        'total': duration_resolution + duration_gen + duration_region,
    }
    for name, value in durations_process.items():
        mpas_ds.attrs['vtx-duration-' + name] = '%.2f' % value

    mpas_ds.to_netcdf(mpas_grid_file)

    if not os.path.isfile(mpas_grid_file):
        raise IOError('Could not update attributes of file ' +
                      mpas_grid_file_tmp)

    os.system('rm -f ' + mpas_grid_file_tmp)

    if do_plots:
        ds = open_mpas_regional_file(mpas_grid_file)

        plot_mpas_darray(ds, 'resolution',
                         cmap='Spectral',
                         ax=None,
                         outfile=path_save + '/resolution_mesh.png',
                         title='<NAME>: <VAR>',
                         name=os.path.basename(mpas_grid_file))

    return mpas_grid_file, graph_info_file
