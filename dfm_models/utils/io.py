"""Code for io and network calls

cody.l.johnson@erdc.dren.mil

"""

from pathlib import Path
from tempfile import NamedTemporaryFile as tmp
from urllib.error import HTTPError
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
import xarray as xr
from pydap.cas.get_cookies import setup_session

from dfm_models._internal import (
    validate_cmems_variable,
    validate_datetime,
    validate_variable,
)


def download_ncoda(lats, lons, t0, tf, variables, region=1, fn=None):
    """Subset HYCOM output using OpenDAP

    :params:
        lats = [south, north] limits of bbox
        lons = [west, east] limits of bbox
        datetime = "Y-M-D HH:mm" string
        variables = list of variables in ["salinity", "water_temp", "surf_el", "water_u", "water_v"]
        region = Hycom re-analysis region (default region=1)

    :returns:
        Xarray Dataset of selected variables

    """

    validate_datetime(t0)
    validate_datetime(tf)
    _variables, surf_el = fix_surf_el(variables)

    # query dataset to get coordinates and convert bbox to indicies for OpenDAP
    coords = xr.open_dataset(request)
    lon_ll = lon2index(lons[0], coords, corr=False)  # lower left longtiude of bbox
    lon_ur = lon2index(lons[1], coords, corr=False)
    lat_ll = lat2index(lats[0], coords)
    lat_ur = lat2index(lats[1], coords)

    request = (
        f"https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ssh?"
        f"lat[lat_ll:1:lat_ur],lon[lon_ll:1:lon_ur],time[0:1:7562],surf_el[0:1:0][lat_ll:1:lat_ur][lon_ll:1:lon_ur]"
    )

    request = request + "".join(
        [
            f"{variable}[0:1:0][0:1:39][{lat_ll}:1:{lat_ur}][{lon_ll}:1:{lon_ur}],"
            for variable in _variables
        ]
    )

    # append surf_el if present
    if surf_el is not None:
        request = (
            request + f"{surf_el}[0:1:0][{lat_ll}:1:{lat_ur}][{lon_ll}:1:{lon_ur}],"
        )

    request = request + "time[0:1:0]"

    ds = xr.open_dataset(request)

    if fn is not None:
        ds.to_netcdf(fn)

    return ds


def download_hycom(lats, lons, datetime, variables, region=1, fn=None):
    """Subset HYCOM output using OpenDAP

    :params:
        lats = [south, north] limits of bbox
        lons = [west, east] limits of bbox
        datetime = "Y-M-D HH:mm" string
        variables = list of variables in ["salinity", "water_temp", "surf_el", "water_u", "water_v"]
        region = Hycom re-analysis region (default region=1)

    :returns:
        Xarray Dataset of selected variables

    """

    datetime = pd.to_datetime(datetime)

    validate_datetime(datetime)

    try:
        validate_variable(variables)
    except NameError:
        raise NameError("Input 'variable' needs to be specified")

    _variables, surf_el = fix_surf_el(variables)

    ymd = datetime.strftime("%Y%m%d")
    hr = datetime.strftime("%H")

    request = (
        f"https://www.ncei.noaa.gov/thredds-coastal/dodsC/hycom_region{region}"
        f"/{ymd}/hycom_glb_regp{region:02d}_{ymd}00_t0{hr}.nc?"
        f"depth[0:1:39],lat[0:1:875],lon[0:1:625]"
    )

    # query dataset to get coordinates and convert bbox to indicies for OpenDAP
    coords = xr.open_dataset(request)
    lon_ll = lon2index(lons[0], coords)  # lower left longtiude of bbox
    lon_ur = lon2index(lons[1], coords)
    lat_ll = lat2index(lats[0], coords)
    lat_ur = lat2index(lats[1], coords)

    request = (
        f"https://www.ncei.noaa.gov/thredds-coastal/dodsC/hycom_region{region}/"
        f"{ymd}/hycom_glb_regp{region:02d}_{ymd}00_t0{hr}.nc?"
        f"depth[0:1:39],lat[{lat_ll}:1:{lat_ur}],lon[{lon_ll}:1:{lon_ur}],"
    )

    request = request + "".join(
        [
            f"{variable}[0:1:0][0:1:39][{lat_ll}:1:{lat_ur}][{lon_ll}:1:{lon_ur}],"
            for variable in _variables
        ]
    )

    # append surf_el if present
    if surf_el is not None:
        request = (
            request + f"{surf_el}[0:1:0][{lat_ll}:1:{lat_ur}][{lon_ll}:1:{lon_ur}],"
        )

    request = request + "time[0:1:0]"

    ds = xr.open_dataset(request)

    if fn is not None:
        ds.to_netcdf(fn)

    return ds


def download_ncom(lats, lons, datetime, variables, region="amseas", fn=None):
    """Subset NCOM output using OpenDAP

    :params:
        lats = [south, north] limits of bbox
        lons = [west, east] limits of bbox
        datetime = "Y-M-D HH:mm" string
        variables = list of variables in ["salinity", "water_temp", "surf_el", "water_u", "water_v"]
        region = Hycom re-analysis region (default region=1)

    :returns:
        Xarray Dataset of selected variables

    """

    datetime = pd.to_datetime(datetime)

    validate_datetime(datetime)

    try:
        validate_variable(variables)
    except NameError:
        raise NameError("Input 'variable' needs to be specified")

    _variables, surf_el = fix_surf_el(variables)

    ymd = datetime.strftime("%Y%m%d")
    hr = datetime.strftime("%H")

    request = (
        f"https://www.ncei.noaa.gov/thredds-coastal/dodsC/{region}/{region}_20130405_to_current/{ymd}"
        f"/ncom_relo_{region}_u_{ymd}00_t0{hr}.nc?lon[0:1:1293],lat[0:1:813],depth[0:1:39]"
    )

    # query dataset to get coordinates and convert bbox to indicies for OpenDAP
    try:
        coords = xr.open_dataset(request)
    except OSError:
        raise OSError(f"URL not found: {request}")

    lon_ll = lon2index(lons[0], coords)  # lower left longtiude of bbox
    lon_ur = lon2index(lons[1], coords)
    lat_ll = lat2index(lats[0], coords)
    lat_ur = lat2index(lats[1], coords)

    request = (
        f"https://www.ncei.noaa.gov/thredds-coastal/dodsC/{region}/{region}_20130405_to_current/{ymd}"
        f"/ncom_relo_{region}_u_{ymd}00_t0{hr}.nc?lon[{lon_ll}:1:{lon_ur}],lat[{lat_ll}:1:{lat_ur}],depth[0:1:39],"
    )

    request = request + "".join(
        [
            f"{variable}[0:1:0][0:1:39][{lat_ll}:1:{lat_ur}][{lon_ll}:1:{lon_ur}],"
            for variable in _variables
        ]
    )

    # append surf_el if present
    if surf_el is not None:
        request = (
            request + f"{surf_el}[0:1:0][{lat_ll}:1:{lat_ur}][{lon_ll}:1:{lon_ur}],"
        )

    request = request + "time[0:1:0]"

    try:
        ds = xr.open_dataset(request)
    except OSError:
        raise OSError(f"URL not found: {request}")

    if fn is not None:
        ds.to_netcdf(fn)

    return ds


def download_ocean_ts(lats, lons, t0, tf, variables, download):
    """Subset NCOM output using OpenDAP

    :params:
        lats = [south, north] limits of bbox
        lons = [west, east] limits of bbox
        datetime = "Y-M-D HH:mm" string
        variables = list of variables in ["salinity", "water_temp", "surf_el", "water_u", "water_v"]
        region = Hycom re-analysis region (default region=1)

    :returns:
        Xarray Dataset of selected variables
    """

    datasets = []

    tmpDir = Path("/tmp") / "".join(
        random.choices(string.ascii_letters + string.digits, k=10)
    )

    tmpDir.mkdir()

    t1 = pd.to_datetime(t0)
    t2 = pd.to_datetime(tf)
    date_range = pd.date_range(t1, t2, freq="3H")

    for datetime in date_range:

        fn = tmpDir / datetime.strftime("%Y%m%d.%H%M.nc")

        try:
            ds = download(lats, lons, datetime, variables, fn=fn)
            datasets.append(ds)

        except OSError:
            print(f"Skipping timestep {datetime} because URL was not found.")
            continue

    print(f"Successfully downloaded data set to: {tmpDir}")

    return xr.concat(datasets, "time")


def download_COOPs(product, station_name, station_id, datum, begin_date, end_date):
    """download predicted water level and append to attribute

    :function: TODO
    :returns: TODO

    """
    request = (
        f"https://tidesandcurrents.noaa.gov/api/datagetter?"
        f"begin_date={begin_date}"
        f"&end_date={end_date}"
        f"&station={station_id}"
        f"&product={product}"
        f"&datum={datum}"
        f"&units=metric"
        f"&time_zone=gmt"
        f"&application=ERDC"
        f"&format=csv"
        f"&interval=h"
    )

    fn = tmp()

    try:
        response, http = urlretrieve(request, fn.name)
    except HTTPError:
        print(f"{station_name} with {station_id} was not found in CO-OPs database.")
        print(f"Check url for errors: {request}")
        raise

    data = pd.read_csv(
        fn,
        index_col=[0],
        parse_dates=True,
        names=["time", "water_level"],
        header=0,
        usecols=[0, 1],
        squeeze=True,
    )

    data.name = f"water_level_{product}"

    return data


def download_nwis(
    station_name, station_id, begin_date, end_date, data_code=60, skiprows=28
):
    """download data from https://nwis.waterdata.usge and outputs as dataframe

    inputs:
    site_name = user specified name for site
    site_no = USGS site number code
    begin_date = first day in timeseries (YYYY-MM-dd)
    end_date = last day in timeseries (YYYY-MM-dd)
    skiprows = number of header rows to skip (default=28)

    return = time series (pandas Series)
    """

    request = (
        f"https://nwis.waterdata.usgs.gov/usa/nwis/uv/"
        f"?cb_{data_code:05d}=on"
        f"&format=rdb&"
        f"site_no={station_id}"
        f"&period="
        f"&begin_date={begin_date}"
        f"&end_date={end_date}"
    )

    fn = tmp()

    try:
        response, http = urlretrieve(request, fn.name)
    except HTTPError:
        print(f"{station_name} with {station_id} was not found in CO-OPs database.")
        print(f"Check url for errors: {request}")
        raise

    # Pandas
    data = pd.read_csv(
        fn,
        sep="\s+",  # noqa: W605
        skiprows=skiprows,
        usecols=[2, 3, 5],
        parse_dates={"datetime_CST": [0, 1]},
        header=0,
        index_col=0,
        names=["date", "time", "data"],
        dtype={"data": float},
        squeeze=True,
    )

    try:
        data.index = (
            data.index.tz_localize("America/Chicago", ambiguous=True)
            .tz_convert("UTC")
            .tz_localize(None)
        )
        data.index = data.index.rename("datetime_UTC")
    except AttributeError as e:  # noqa: F841
        print("Problem converting datetime to UTC. Check data")
        return np.nan

    return data


def download_cmems_ts(lats, lons, t0, tf, variables, fn=None):
    """Subset CMEMS output using OpenDAP

    :params:
        lats = [south, north] limits of bbox
        lons = [west, east] limits of bbox
        t0 = datetime for start of time series
        tf = datetime for end of time series
        variables = list of variables in ["zos", "uo", "vo", "so", "thetao"]

    :returns:
        Xarray Dataset of selected variables
    """

    validate_datetime(t0)
    validate_datetime(tf)

    try:
        validate_cmems_variable(variables)
    except NameError:
        raise NameError("Input 'variable' needs to be specified")

    _variables, zos = fix_zos(variables)

    request = (
        "https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m?"
        "longitude[0:1:4319],latitude[0:1:2040],depth[0:1:49],time[0:1:10012]"
    )

    # query dataset to get coordinates and convert bbox to indicies for OpenDAP
    coords = xr.open_dataset(request)
    lon_ll = cmemslon2index(lons[0], coords)  # lower left longtiude of bbox
    lon_ur = cmemslon2index(lons[1], coords)
    lat_ll = cmemslat2index(lats[0], coords)
    lat_ur = cmemslat2index(lats[1], coords)
    t0i = time2index(t0, coords)
    tfi = time2index(tf, coords)

    request = (
        f"https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m?"
        f"longitude[{lon_ll}:1:{lon_ur}],latitude[{lat_ll}:1:{lat_ur}],depth[0:1:49],time[{t0i}:1:{tfi}],"
    )

    request = request + "".join(
        [
            f"{variable}[{t0i}:1:{tfi}][0:1:49][{lat_ll}:1:{lat_ur}][{lon_ll}:1:{lon_ur}]"
            for variable in _variables
        ]
    )

    # append surf_el if present
    if zos is not None:
        request = (
            request + f"{zos}[{t0i}:1:{tfi}][{lat_ll}:1:{lat_ur}][{lon_ll}:1:{lon_ur}]"
        )

    ds = xr.open_dataset(request)

    if fn is not None:
        ds.to_netcdf(fn)

    return ds


def copernicusmarine_session(username, password):
    cas_url = "https://cmems-cas.cls.fr/cas/login"
    session = setup_session(cas_url, username, password)
    session.cookies.set("CASTGC", session.cookies.get_dict()["CASTGC"])


###########
# helpers #
###########
def fix_surf_el(variables):
    """change variables if surf_el is contained"""

    if "surf_el" in set(variables):

        _variables = variables.copy()

        _variables.remove("surf_el")

        return _variables, "surf_el"

    else:

        return variables, None


def fix_zos(variables):
    """change variables if surf_el is contained"""

    if "zos" in set(variables):

        _variables = variables.copy()

        _variables.remove("zos")

        return _variables, "zos"

    else:

        return variables, None


def time2index(t, coords):
    """convert time to index for OpenDAP request"""
    validate_datetime(t)
    time = pd.to_datetime(coords.time.values)
    return np.argmin(np.abs(t - time))


def lon2index(lon, coords, corr=True):
    """convert longitude to index for OpenDAP request"""
    if corr:
        if lon < 0:
            lon += 360
    lons = coords.lon.values
    return np.argmin(np.abs(lons - lon))


def lat2index(lat, coords):
    """convert latitude to index for OpenDAP request"""

    lats = coords.lat.values

    return np.argmin(np.abs(lats - lat))


def cmemslon2index(lon, coords):
    """convert longitude to index for OpenDAP request"""
    lons = coords.longitude.values
    return np.argmin(np.abs(lons - lon))


def cmemslat2index(lat, coords):
    """convert latitude to index for OpenDAP request"""

    lats = coords.latitude.values
    return np.argmin(np.abs(lats - lat))
