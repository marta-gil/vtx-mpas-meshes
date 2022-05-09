import argparse
import os
import xarray as xr

from vtxmpasmeshes.plot_utilities import view_mpas_regional_mesh, \
    plot_expected_resolution_rings

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument(
    "-g", "--grid", type=str, required=True,
    help="Name of an MPAS grid.nc",
)

parser.add_argument(
    "-o", "--outfile", type=str, default=None,
    help="File to save the MPAS plot",
)

args = parser.parse_args()

if not os.path.exists(args.grid):
    raise IOError('File does not exist: ' + args.grid)

view_mpas_regional_mesh(args.grid, outfile=args.outfile)
