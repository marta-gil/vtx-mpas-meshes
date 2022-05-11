import xarray as xr
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from geopy.distance import distance


# OPEN AND CLOSE PLOTS

def start_cartopy_map_axis(zorder=1):
    ax = plt.axes(projection=ccrs.PlateCarree())  # projection type
    add_cartopy_details(ax, zorder=zorder)
    return ax


def add_cartopy_details(ax, zorder=1):
    ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=zorder)
    ax.coastlines(resolution='10m', zorder=zorder)

    gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle='--',
                      zorder=zorder)
    gl.top_labels = False
    gl.right_labels = False


def close_plot(fig=None, size_fig=None, pdf=None, outfile=None,
               force_show=False):

    if size_fig is None:
        size_fig = [10, 8]

    if fig is None:
        fig = plt.gcf()
    fig.set_size_inches(size_fig)

    if outfile is not None:
        plt.savefig(outfile)

    if pdf is not None:
        pdf.savefig(fig)

    if (outfile is None and pdf is None) or force_show:
        plt.show()

    plt.close()


# HANDLE PLOT KWARGS FROM ALL PASSED KWARGS

def set_plot_kwargs(da=None, list_darrays=None, **kwargs):
    plot_kwargs = {k: v for k, v in kwargs.items()
                   if k in ['cmap', 'vmin', 'vmax']
                   and v is not None}

    if 'cmap' not in plot_kwargs:
        plot_kwargs['cmap'] = 'Spectral'

    vmin = plot_kwargs.get('vmin', None)
    if vmin is None:
        if da is not None:
            vmin = np.min(da)
        elif list_darrays is not None:
            vmin = np.min([v.min() for v in list_darrays if v is not None])
    if vmin is not None:
        plot_kwargs['vmin'] = vmin

    vmax = plot_kwargs.get('vmax', None)
    if vmax is None:
        if da is not None:
            vmax = np.max(da)
        elif list_darrays is not None:
            vmax = np.max([v.max() for v in list_darrays if v is not None])

    if vmax is not None:
        plot_kwargs['vmax'] = vmax

    return plot_kwargs


# BORDERS OF THE PLOT

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


def get_borders_at_distance(distance_km, centerlat=0., centerlon=0.):
    len_grid = distance(kilometers=distance_km)

    maxlat = len_grid.destination(point=(centerlat, centerlon),
                                  bearing=0).latitude
    minlat = len_grid.destination(point=(centerlat, centerlon),
                                  bearing=180).latitude
    maxlon = len_grid.destination(point=(centerlat, centerlon),
                                  bearing=90).longitude
    minlon = len_grid.destination(point=(centerlat, centerlon),
                                  bearing=270).longitude

    return minlon, maxlon, minlat, maxlat


# BASIC LATLON PLOT

def plot_latlon_cartopy(darray, ax=None, title='', borders=None, **kwargs):
    # For a DataArray that has lat/lon coordinates

    # Start the axis
    if ax is None:
        ax = start_cartopy_map_axis(zorder=4)

    # Set the extent
    if borders is None:
        borders = find_borders(darray.lat.values, darray.lon.values)
        # borders = [minlon, maxlon, minlat, maxlat]
    ax.set_extent(borders, crs=ccrs.PlateCarree())

    # Plot the array (automatic xarray plot)
    darray.plot(ax=ax, **kwargs)  # plot the resolution

    ax.set_title(title)
    return


# MPAS PLOT

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


def plot_mpas_darray(ds, vname, ax=None, outfile=None, **kwargs):

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
        ax = start_cartopy_map_axis()

    if 'borders' not in kwargs:
        lats = ds['latitude'].values.flatten()
        lons = ds['longitude'].values.flatten()
        borders = find_borders(lats, lons)
    else:
        borders = kwargs['borders']

    ax.set_extent(borders, crs=ccrs.PlateCarree())

    plot_kwargs = set_plot_kwargs(da=da, **kwargs)

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

        color = colorvalue(value, **plot_kwargs)

        ax.fill(lons, lats, edgecolor=None, linewidth=0.0,
                facecolor=color)

    units = da.attrs.get('units', '')
    name = kwargs.get('name', '')
    ncells = str(len(da.values.flatten()))
    title = kwargs.get('title', '')
    title = title.replace('<VAR>', vname).replace('<UNITS>', units)
    title = title.replace('<NAME>', name).replace('<NCELLS>', ncells)
    ax.set_title(title)

    if final:
        title_legend = kwargs.get('title_legend', '<VAR>: <UNITS>')
        title_legend = title_legend.replace('<VAR>', vname)
        title_legend = title_legend.replace('<UNITS>', units)
        add_colorbar(ax, label=title_legend, **plot_kwargs)

        close_plot(outfile=outfile)
    return


# TOOLS FOR MULTIPLE PLOTS (Also can be used for 1 plot)

def add_colorbar(axs, fig=None, label=None, **plot_kwargs):
    if fig is None:
        fig = plt.gcf()

    try:
        x = axs[0, 0]
    except:
        try:
            x = axs[0]
            n = len(axs)
        except:
            axs = np.array([axs]).reshape([1, 1])
        else:
            axs = axs.reshape([n, 1])

    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(vmin=plot_kwargs['vmin'],
                                      vmax=plot_kwargs['vmax'], clip=True),
            cmap=plot_kwargs['cmap']),
        ax=axs[:, :], shrink=0.6)
    cbar.ax.locator_params(nbins=10)
    if label is not None:
        cbar.set_label(label)

    return


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


def get_max_borders(list_results, margin='factor2', namelat='latitude',
                    namelon='longitude'):
    borders = np.array([find_borders(lats=ds[namelat].values,
                                     lons=ds[namelon].values,
                                     margin=margin)
                        for ds in list_results])
    max_borders = borders[:, 0].min(), borders[:, 1].max(), \
                  borders[:, 2].min(), borders[:, 3].max()
    return np.array(max_borders)
