import numpy as np
import xarray as xr
import math


derived_variables = {
    'latCell': ['latitude'],
    'lonCell': ['longitude'],
    'latVertex': ['latitudeVertex'],
    'lonVertex': ['longitudeVertex'],
    'areaCell': ['area', 'resolution'],
}


def distance_latlon_matrix(lat, lon, lat_ref=0., lon_ref=0., do_tile=False):
    if do_tile:
        lat = np.tile(lat, (len(lon), 1)).transpose()
        lon = np.tile(lon, (len(lat), 1))

    lon, lat, lon_ref, lat_ref = map(np.radians, [lon, lat, lon_ref, lat_ref])

    newlon = lon_ref - lon
    newlat = lat_ref - lat

    haver_formula = np.sin(newlat / 2.0) ** 2 + np.cos(lat) * \
                    np.cos(lat_ref) * np.sin(newlon / 2.0) ** 2

    dist = 2 * np.arcsin(np.sqrt(haver_formula))
    km = 6367 * dist  # 6367 for distance in KM for miles use 3958
    return km


def add_mpas_mesh_variables(ds, full=True):
    for v in ds.data_vars:
        if v not in derived_variables:
            continue

        newvs = derived_variables[v]

        for newv in newvs:
            if newv in ds:
                #print(newv + ' already here')
                continue

            if 'lat' in v or 'lon' in v:
                ds[newv] = xr.apply_ufunc(np.rad2deg, ds[v])
                ds[newv] = ds[newv].where(ds[newv] <= 180.0, ds[newv] - 360.0)
                ds[newv].attrs['units'] = 'degrees'

            elif newv == 'area':
                radius_circle = ds.attrs.get('sphere_radius', 1.0)
                if radius_circle == 1:
                    #print('need to correct to earth radius!!')
                    correction_rad_earth = 6371220.0
                else:
                    correction_rad_earth = 1

                ds[newv] = (ds[v] / 10 ** 6) * correction_rad_earth**2
                ds[newv].attrs['units'] = 'km^2 (assuming areaCell in m^2)'
                ds[newv].attrs['long_name'] = 'Area of the cell in km^2'

            elif newv == 'resolution':
                radius_circle = ds.attrs.get('sphere_radius', 1.0)
                if radius_circle == 1:
                    #print('need to correct to earth radius!!')
                    correction_rad_earth = 6371220.0
                else:
                    correction_rad_earth = 1

                # km^2 (assuming areaCell in m^2)
                area = (ds[v] / 10 ** 6) * correction_rad_earth**2

                ds[newv] = 2 * (xr.apply_ufunc(np.sqrt, area / math.pi))
                ds[newv].attrs['units'] = 'km'
                ds[newv].attrs['long_name'] = 'Resolution of the cell (approx)'

    if full:
        ds = add_distance_to_reference(ds)
        ds = compute_metrics_edge_lengths(ds)

    return ds


def get_center(lats, lons):
    """
    Given 2 DataArrays of lats, lons, it returns the grid center lat,
    lon.
    :param lats: DataArray of latitudes
    :param lons: DataArray of longitudes
    :return: lat, lon: floats
    """

    lat = float(lats[abs(lats - lats.mean()).argmin()])
    lon = float(lons[abs(lons - lons.mean()).argmin()])

    return lat, lon


def get_distance_to_center(lats, lons,
                           center_lat=None, center_lon=None):
    """
    Using the haversine formula (spherical distance), it returns the
    distance in km of each lat, lon point to the
    central point.
    """
    if center_lat is None or center_lon is None:
        center_lat, center_lon = get_center(lats, lons)

    radius = 6371.  # km
    d_lat = np.radians(lats - center_lat)  # lat distance in radians
    d_lon = np.radians(lons - center_lon)  # lon distance in radians

    a = (np.sin(d_lat / 2.) * np.sin(d_lat / 2.) +
         np.cos(np.radians(center_lat)) * np.cos(np.radians(lats)) *
         np.sin(d_lon / 2.) * np.sin(d_lon / 2.))
    c = 2. * np.arctan2(np.sqrt(a), np.sqrt(1. - a))
    d = radius * c

    return d


def compute_metrics_edge_lengths(ds):
    distances = []

    for i, cell in enumerate(ds['nCells'].values):
        vals = ds['verticesOnCell'].sel(nCells=cell).values

        num_sides = int(ds['nEdgesOnCell'].sel(nCells=cell))
        vals = vals[:num_sides] - 1
        lats = ds['latitudeVertex'].values[vals]
        lons = ds['longitudeVertex'].values[vals]
        lat_ref = lats[-1]
        lon_ref = lons[-1]

        my_cell_dists = []
        for i in range(ds.dims['maxEdges']):
            if i >= num_sides:
                my_cell_dists.append(np.nan)
                continue
            d = distance_latlon_matrix(lats[i], lons[i],
                                       lat_ref=lat_ref,
                                       lon_ref=lon_ref, do_tile=False)

            my_cell_dists.append(d)
            lat_ref, lon_ref = lats[i], lons[i]
        distances.append(my_cell_dists)

    ds['edgesLength'] = xr.DataArray(data=distances,
                                     dims=('nCells', 'maxEdges'))
    ds['edgesLength'].attrs = {
        'name': 'Length of edges',
        'units': 'km',
        'long_name': 'Haversine distance between adjacent vertices of a cell.'
    }

    ds['mean_edgesLength'] = ds['edgesLength'].mean(dim='maxEdges')
    ds['mean_edgesLength'].attrs = {
        'name': 'Mean edge length',
        'units': 'km',
        'long_name': 'Mean edge length of a cell.'
    }

    ds['min_edgesLength'] = ds['edgesLength'].min(dim='maxEdges')
    ds['min_edgesLength'].attrs = {
        'name': 'Minimum edge length',
        'units': 'km',
        'long_name': 'Minimum edge length of a cell.'
    }

    ds['max_edgesLength'] = ds['edgesLength'].max(dim='maxEdges')
    ds['max_edgesLength'].attrs = {
        'name': 'Maximum edge',
        'units': 'km',
        'long_name': 'Maximum edge length of a cell.'
    }

    ds['rmse_edgesLength'] = xr.apply_ufunc(np.sqrt,
                                            (ds['edgesLength'] ** 2).mean(
                                                dim='maxEdges'))
    ds['rmse_edgesLength'].attrs = {
        'name': 'Rmse edge length',
        'units': 'km',
        'long_name': 'Root mean squared edge length.'
    }

    ds['ration_edgesLength'] = ds['min_edgesLength'] / ds['max_edgesLength']
    ds['ration_edgesLength'].attrs = {
        'name': 'Ratio of edge lengths',
        'units': '',
        'long_name': 'Ratio between minimum and maximum edge lengths '
                     'of a cell.'
    }

    ds['cellDistortion'] = xr.apply_ufunc(np.sqrt, (
                (ds['edgesLength'] - ds['rmse_edgesLength']) ** 2).mean(
        dim='maxEdges')) / ds['rmse_edgesLength']
    ds['cellDistortion'].attrs = {
        'name': 'Cell Distortion',
        'units': '',
        'long_name': 'Cell Distortion: Haversine distance between adjacent '
                     'vertices of a cell. '
                     'https://pedrosp.ime.usp.br/papers/PeixotoBarros2013.pdf'
    }
    return ds


def add_distance_to_reference(ds, **kwargs):
    lat_ref = kwargs.get('lat_ref', ds.attrs['vtx-param-lat_ref'])
    lon_ref = kwargs.get('lon_ref', ds.attrs['vtx-param-lon_ref'])

    d = get_distance_to_center(ds['latitude'].values, ds['longitude'].values,
                               center_lat=lat_ref, center_lon=lon_ref)

    ds['cellDistance'] = xr.DataArray(data=d, dims=('nCells'))
    ds['cellDistance'].attrs = {
        'name': 'Distance to ref. point',
        'units': 'km',
        'long_name': 'Distance from cell center to reference point'
    }
    return ds


def open_mpas_regional_file(file, **kwargs):
    ds = xr.open_dataset(file)
    ds = add_mpas_mesh_variables(ds, **kwargs)
    return ds

