import argparse
import os

from vtxmpasmeshes.mpas_plots import compare_plot_mpas_regional_meshes

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-gs', '--grids', default=[],
                    help='names of the grids',
                    type=lambda s: [item for item in s.split(',')]
                    )

parser.add_argument(
    "-o", "--outfile", type=str, default=None,
    help="File to save the MPAS plot",
)

args = parser.parse_args()

DATA_FOLDER = '/home/marta/PycharmProjects/vtx-mpas-meshes/data'

grids = []
for g in args.grids:
    grid = DATA_FOLDER + '/' + g + '/' + g + '.grid.nc'
    if not os.path.exists(grid):
        raise IOError('File does not exist: ' + grid)
    grids.append(grid)

compare_plot_mpas_regional_meshes(grids,
                                  outfile=args.outfile,
                                  suptitle='Meshes comparison',
                                  each_title='<NAME>: <NCELLS> cells')
