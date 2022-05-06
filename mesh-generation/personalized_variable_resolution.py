import argparse

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def distance_latlon_matrix(lat, lon, lat_ref=0., lon_ref=0., do_tile=False):
    if do_tile:
        lat = np.tile(lat, (len(lon), 1)).transpose()
        lon = np.tile(lon, (len(lat), 1))

    lon, lat, lon_ref, lat_ref = map(np.radians, [lon, lat, lon_ref, lat_ref])

    newlon = lon_ref - lon
    newlat = lat_ref - lat

    haver_formula = np.sin(newlat / 2.0) ** 2 + np.cos(lat) * np.cos(lat_ref) * np.sin(newlon / 2.0) ** 2

    dist = 2 * np.arcsin(np.sqrt(haver_formula))
    km = 6367 * dist  # 6367 for distance in KM for miles use 3958
    return km


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


def find_borders(lats, lons, margin='factor2'):
    lats = lats.flatten()
    lons = lons.flatten()

    limits = lons.min(), lons.max(), lats.min(), lats.max()
    if 'factor' in margin:
        delta = int(margin[-1])
        deltalat = np.abs(limits[1] - limits[0]) / (10 * delta)
        deltalon = np.abs(limits[3] - limits[2]) / (10 * delta)
        minlon = limits[0] - deltalon
        maxlon = limits[1] + deltalon
        minlat = limits[2] - deltalat
        maxlat = limits[3] + deltalat
        limits = minlon, maxlon, minlat, maxlat
    elif 'cells' in margin:
        delta = int(margin[-1])
        deltalat = np.abs(limits[1] - limits[0])
        deltalon = np.abs(limits[3] - limits[2])
        minlon = limits[0] - deltalon * delta
        maxlon = limits[1] + deltalon * delta
        minlat = limits[2] - deltalat * delta
        maxlat = limits[3] + deltalat * delta
        limits = minlon, maxlon, minlat, maxlat

    return limits


def plot_cartopy(darray, ax=None, title='',  borders=None, **kwargs):
    if ax is None:
        ax = plt.axes(projection=ccrs.PlateCarree())

    if borders is None:
        borders = find_borders(darray.lat.values, darray.lon.values)
        print(borders)

    ax.set_extent(borders, crs=ccrs.PlateCarree())

    darray.plot(ax=ax, **kwargs)  # plot the resolution

    ax.set_title(title)

    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.coastlines(resolution='10m')
    gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    return


def view_resolution_map(ds, pdfname=None, list_distances=None):

    kwargs = {'cmap': 'Spectral', 'vmin': 0, 'levels': 41}

    if list_distances is None:
        list_distances = [1000, 500, 200, 50]

    pdf = None
    if pdfname is not None:
        pdf = PdfPages(pdfname)

    # Region
    # we don't want to plot too much
    region = ds.where(ds['distance'] <= max(list_distances))
    region = region.dropna('lat', how='all').dropna('lon', how='all')
    my_zero = region['distance'].argmin(['lat', 'lon'])
    mylat = int(my_zero['lat'])
    mylon = int(my_zero['lon'])

    axis = region.isel(lat=mylat, lon=range(mylon, region.dims['lon']))
    axis = axis.squeeze(drop=True)
    axis = axis.assign_coords({'distance': axis['distance']})
    axis = axis.swap_dims({'lon': 'distance'})

    for di in list_distances:
        print('\t .. plotting for distances <= %.0fkm' % di)
        plt.subplot(121)
        axis['resolution'].where(axis['distance'] <= di).plot()
        plt.title('Radial Resolution')

        x = region['resolution'].where(region['distance'] <= di)
        x = x.dropna('lat', how='all').dropna('lon', how='all')
        ax = plt.subplot(122, projection=ccrs.PlateCarree())
        plot_cartopy(x, ax=ax, title='Resolution map', **kwargs)

        fig = plt.gcf()
        fig.suptitle('Resolution (km). Distance closer than %.0fkm' % di)
        if pdf is not None:
            pdf.savefig(fig)
        else:
            plt.show()
        plt.close()

    if pdf is not None:
        pdf.close()

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

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
                     again until reaching 100km (to save space).
                     The requested MPAS region should be circular and 
                     have a radius of <size>+<margin>. The buffer 
                     generated by the MPAS-Limited-Area code will then
                     consist of a few "rings" of <lowresolution> cells.
                     \n
            """
    )

    parser.add_argument(
        "-highr", "--highresolution", required=True, type=float,
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
        "-clat", "--clat", default=0., type=float,
        help="Central latitude.",
    )

    parser.add_argument(
        "-clon", "--clon", default=0., type=float,
        help="Central longitude.",
    )

    # -p generates plots
    parser.add_argument(
        "-p", "--withplots", action="store_true",
        help="generate plots to view the resolution.",
    )

    parser.add_argument(
        "-pdf", "--pdfname", default=None, type=str,
        help="pdffile where to save the resolution plots.",
    )

    args = parser.parse_args()

    ds = variable_resolution_latlonmap(args.grid,
                                       highresolution=args.highresolution,
                                       lowresolution=args.lowresolution,
                                       size=args.size,
                                       margin=args.margin,
                                       lat_ref=args.clat,
                                       lon_ref=args.clon,
                                       )

    print(ds)

    if args.withplots or args.pdfname is not None:
        print('Plotting')
        view_resolution_map(ds, pdfname=args.pdfname,
                            list_distances=[1000, 500,
                                            ds.attrs['border'],
                                            ds.attrs['radius'],
                                            ])

