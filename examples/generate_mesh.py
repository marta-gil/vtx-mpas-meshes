import argparse
import os

from vtxmpasmeshes.mesh_generator import full_generation_process

DATA_FOLDER = '/home/marta/PycharmProjects/vtx-mpas-meshes/data'

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