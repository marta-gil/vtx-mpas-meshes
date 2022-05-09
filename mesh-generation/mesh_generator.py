import argparse
import os
import time

import xarray as xr

from personalized_variable_resolution import variable_resolution_latlonmap, \
    view_resolution_map

from jigsaw_generator import jigsaw_gen_sph_grid

from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.mesh.conversion import convert
from mpas_tools.io import write_netcdf


PATH_LIMITED_AREA = '/home/marta/PycharmProjects/MPAS-Limited-Area'
DATA_FOLDER ='/home/marta/PycharmProjects/vtx-mpas-meshes/data'


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
        path_save = os.path.dirname(mpas_grid_file)
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
        mpas_ds.attrs['vtx-param-' + name] = value

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

    return mpas_grid_file, graph_info_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "-n", "--name", default="doughnut", type=str,
        help="output basename for directory and files."
    )

    parser.add_argument(
        "-g", "--grid", type=str, default='doughnut',
        help="""
            Grid option: \n 
            "  doughnut": 
                     <High resolution> area of a certain radius <size>. 
                     Linear increase of resolution to a certain <low
                     resolution> value after <margin>km. 
                     Keep the constant low resolution value for a while
                     (10*low_resolution)km and then increase it linearly
                     again until reaching 1000km (to save space).
                     The requested MPAS region should be circular and 
                     have a radius of <size>+<margin>. The buffer 
                     generated by the MPAS-Limited-Area code will then
                     consist of a few "rings" of <lowresolution> cells.
                     \n
            """
    )

    parser.add_argument(
        "-highr", "--highresolution", default=3, type=float,
        help="Highest-resolution of grid (km).",
    )

    parser.add_argument(
        "-lowr", "--lowresolution", default=None, type=float,
        help="Lowest-resolution of grid (km).",
    )

    parser.add_argument(
        "-size", "--size", default=None, type=float,
        help="Radius of the highest-resolution area of the grid (km).",
    )

    parser.add_argument(
        "-margin", "--margin", default=None, type=float,
        help="Size of the variable resolution boundary around the "
             "high resolution area (km).",
    )

    parser.add_argument(
        "-lat", "--lat_ref", default=0., type=float,
        help="Central latitude.",
    )

    parser.add_argument(
        "-lon", "--lon_ref", default=0., type=float,
        help="Central longitude.",
    )

    # -p generates plots
    parser.add_argument(
        "-p", "--withplots", action="store_true",
        help="generate plots to view the mesh.",
    )

    # -o overwrite
    parser.add_argument(
        "-o", "--overwrite", action="store_true",
        help="overwrite existing folder.",
    )

    args = parser.parse_args()

    if args.name == '':
        raise ValueError('Please give a non trivial name.')

    folder = DATA_FOLDER + '/' + args.name + '/'

    if os.path.isdir(folder):
        if not args.overwrite:
            print('Sure to overwrite?')
            raise IOError('For security, overwriting is disabled. Give '
                          'different tests different names or erase the'
                          'existing folder: ' + folder)
        else:
            print('Overwriting folder ' + folder)
            os.system('rm -rf ' + folder)

    os.system('mkdir -p ' + folder)
    basename = folder + args.name

    mesh_file, mesh_graph_info = full_generation_process(
        basename + '.grid.nc',
        args.grid,
        redo=args.overwrite,
        do_plots=args.withplots,
        highresolution=args.highresolution,
        lowresolution=args.lowresolution,
        size=args.size,
        margin=args.margin,
        lat_ref=args.lat_ref,
        lon_ref=args.lon_ref,
    )

    print('\n' + '*' * 30)
    print('\nDONE. This is the mesh ' + mesh_file)
