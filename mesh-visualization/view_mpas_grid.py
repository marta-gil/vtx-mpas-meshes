import argparse
import os
import xarray as xr


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "-g", "--grid", type=str, required=True,
        help="Name of an MPAS grid.nc",
    )

    args = parser.parse_args()

    if not os.path.exists(args.grid):
        raise IOError('File does not exist: ' + args.grid)

    mpas_ds = xr.open_dataset(args.grid)
    print(mpas_ds)

    # TODO - MPAS plots
