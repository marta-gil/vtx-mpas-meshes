import os.path

import xarray as xr
import numpy as np
import shapely

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.geodesic import Geodesic

from vtxmpasmeshes.dataset_utilities import open_mpas_regional_file


mpas_mesh_src_plots = ['latitude', 'longitude', 'latitudeVertex',
                       'longitudeVertex', 'verticesOnCell']


def are_mpas_cells_plottable(ds):
    if 'nCells' not in ds.dims:
        return False

    varsds = list(ds.data_vars.keys())
    if all(x in varsds for x in mpas_mesh_src_plots):
        return True
    else:
        return False


def colorvalue(val, cmap='Spectral', vmin=None, vmax=None):
    """
    Given a value and the range max, min, it returns the associated
    color of the desired cmap.
    :param val: float
    :param cmap: str
    :param vmin: float (default None)
    :param vmax: float (default None)
    :return: cm(norm_val): color
    """
    # Get a colormap instance, defaulting to rc values if name is None.
    cm = mpl.cm.get_cmap(cmap, None)
    if vmin is None:
        vmin = xr.DataArray.min().values  # min value of the array
    if vmax is None:
        vmax = xr.DataArray.max().values  # max value of the array
    if vmin == vmax:
        # A class which, when called, linearly normalizes data into the
        # [0.0, 1.0] interval.
        norm_val = mpl.colors.Normalize(vmin=vmin - 1, vmax=vmax + 1,
                                        clip=True)(val)
    else:
        norm_val = mpl.colors.Normalize(vmin=vmin, vmax=vmax,
                                        clip=True)(val)
    return cm(norm_val)


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


def get_max_borders(list_results, margin='factor2'):
    borders = np.array([r.get_borders(margin=margin) for r in list_results])
    max_borders = borders[:, 0].min(), borders[:, 1].max(), \
                  borders[:, 2].min(), borders[:, 3].max()
    return np.array(max_borders)


def add_colorbar(axs, fig, label=None, **plot_kwargs):
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(vmin=plot_kwargs['vmin'],
                                      vmax=plot_kwargs['vmax'], clip=True),
            cmap=plot_kwargs['cmap']),
        ax=axs[:, :], shrink=0.6)
    cbar.ax.locator_params(nbins=10)
    if label is not None:
        cbar.set_label(label)


def get_plot_size(numplots, nrows=None, ncols=None):
    if nrows is not None and ncols is not None:
        return nrows, ncols

    if nrows is not None:
        ncols = int((numplots + 1) // nrows)
    elif ncols is not None:
        nrows = int((numplots + 1) // ncols)
    else:
        if numplots <= 3:
            nrows, ncols = 1, numplots
        elif numplots == 4:
            nrows, ncols = 2, 2
        elif numplots <= 10:
            ncols = 3
            nrows = int((numplots + 1) // ncols)
        else:
            ncols = 4
            nrows = int((numplots + 1) // ncols)

    return nrows, ncols


def plot_cartopy(darray, ax=None, title='',  borders=None, **kwargs):
    if ax is None:
        ax = plt.axes(projection=ccrs.PlateCarree())

    if borders is None:
        borders = find_borders(darray.lat.values, darray.lon.values)

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


def set_plot_kwargs(da=None, list_darrays=None, **kwargs):
    plot_kwargs = {k: v for k, v in kwargs.items()
                   if k in ['cmap', 'title', 'vmin', 'vmax']}

    if 'title' not in plot_kwargs:
        plot_kwargs['title'] = '<NAME>: <VAR>'

    if 'cmap' not in plot_kwargs:
        plot_kwargs['cmap'] = 'Spectral'

    vmin, vmax = None, None
    if da is not None:
        vmin = np.min(da)
        vmax = np.max(da)
    elif list_darrays is not None:
        vmin = np.min([v.min() for v in list_darrays if v is not None])
        vmax = np.max([v.max() for v in list_darrays if v is not None])

    if 'vmin' not in plot_kwargs:
        if vmin is not None:
            plot_kwargs['vmin'] = vmin

    if 'vmax' not in plot_kwargs:
        if vmax is not None:
            plot_kwargs['vmax'] = vmax
    return plot_kwargs


def view_mpas_regional_mesh(mpas_grid_file, outfile=None, **kwargs):

    ds = open_mpas_regional_file(mpas_grid_file)

    print(ds)

    # PLOTS
    plot_kwargs = set_plot_kwargs(da=ds['resolution'], **kwargs)
    ax = plt.axes(projection=ccrs.PlateCarree())
    axs = np.array([ax]).reshape([1, 1])

    plot_mpas_darray(ds, 'resolution', ax=axs[0, 0], **plot_kwargs,
                     name=os.path.basename(mpas_grid_file))

    add_colorbar(axs, plt.gcf(), label='Resolution (km)',
                 **plot_kwargs)

    plot_expected_resolution_rings(ds, ax=axs[0, 0])

    if outfile is not None:
        plt.savefig(outfile)
    else:
        plt.show()
    plt.close()

    return ds


def plot_expected_resolution_rings(ds, rings=None, outfile=None, ax=None):
    final = False
    if ax is None:
        final = True
        ax = plt.axes(projection=ccrs.PlateCarree())  # projection type
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.coastlines(resolution='10m')

        gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

    lat = ds.attrs['vtx-param-lat_ref']
    lon = ds.attrs['vtx-param-lon_ref']

    if rings is None:
        rings = ['size', 'radius', 'border']

    gd = Geodesic()
    for ring in rings:
        rad = ds.attrs['vtx-param-' + ring]
        cp = gd.circle(lon=lon, lat=lat, radius=rad * 1000)
        geom = shapely.geometry.Polygon(cp)
        ax.add_geometries((geom,), crs=ccrs.PlateCarree(),
                          facecolor='none', edgecolor='black',
                          linewidth=1.5)

    if final:
        fig = plt.gcf()
        fig.set_size_inches([10, 8])

        if outfile is not None:
            plt.savefig(outfile)

        plt.show()
        plt.close()

    return


def plot_mpas_darray(ds, vname, cmap='Spectral', ax=None, outfile=None,
                     title='', **kwargs):

    if vname not in ds.data_vars:
        print('Unplottable Data Array ' + vname)
        print(ds)
        return

    da = ds[vname]
    for coord in ['time', 'lev']:
        if coord in da.dims:
            print('Selecting first slice for ' + coord + '.')
            da = da.isel({coord: 0})

    final = False
    if ax is None:
        final = True
        ax = plt.axes(projection=ccrs.PlateCarree())  # projection type

    if 'borders' not in kwargs:
        lats = ds['latitude'].values.flatten()
        lons = ds['longitude'].values.flatten()
        borders = find_borders(lats, lons)
    else:
        borders = kwargs['borders']

    ax.set_extent(borders, crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.coastlines(resolution='10m')

    gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    min_val = kwargs.get('vmin', da.min().values)
    max_val = kwargs.get('vmax', da.max().values)

    numvert = float(ds.dims['nVertices'])
    last_val = float(ds['verticesOnCell'].sel(maxEdges=-1).values.mean())

    if abs(last_val - numvert) < last_val / 2:
        fillval = numvert
    else:
        fillval = 0.0

    for i, cell in enumerate(da['nCells'].values):
        value = da.sel(nCells=cell)

        vals = ds['verticesOnCell'].sel(nCells=cell).values

        vals = vals[vals != fillval] - 1

        lats = ds['latitudeVertex'].sel(nVertices=vals)
        lons = ds['longitudeVertex'].sel(nVertices=vals)

        color = colorvalue(value, cmap=cmap, vmin=min_val,
                           vmax=max_val)

        ax.fill(lons, lats, edgecolor=None, linewidth=0.0,
                facecolor=color)

    units = da.attrs.get('units', '')
    title = title.replace('<VAR>', vname).replace('<UNITS>', units)
    name = kwargs.get('name', '')
    ncells = str(len(da.values.flatten()))
    title = title.replace('<NAME>', name).replace('<NCELLS>', ncells)
    ax.set_title(title)

    if final:
        fig = plt.gcf()
        fig.set_size_inches([10, 8])

        cbar = fig.colorbar(  # get colorbar
            mpl.cm.ScalarMappable(
                norm=mpl.colors.Normalize(vmin=min_val, vmax=max_val,
                                          clip=True),
                cmap=cmap),
            ax=ax)
        cbar.ax.locator_params(nbins=10)
        title_legend = '<VAR> (<UNITS>)'
        title_legend = title_legend.replace('<VAR>',
                                            vname).replace('<UNITS>', units)
        cbar.set_label(title_legend)

        if outfile is not None:
            plt.savefig(outfile)

        plt.show()
        plt.close()

    return

