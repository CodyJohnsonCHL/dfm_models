"""Object representing FM model

cody.l.johnson@erdc.dren.mil

"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator


def interp2mesh(da, mesh_data, mod180=False):
    """Interpolate variabel from Xarray DataArray to grid

    :params:
        da = DataArray of variable to interpolate
        mesh_data = Xarray DataSet representing UGRID

    :returns: TODO

    """
    xv = mesh_data.NetNode_x.values
    xv[xv < 0] = xv[xv < 0] + 360
    yv = mesh_data.NetNode_y.values
    pv = np.column_stack((xv, yv))

    x = da.lon.values

    if mod180:
        x[x < 0] = x[x < 0] + 360

    y = da.lat.values
    X, Y = np.meshgrid(x, y)
    x = X.ravel()
    y = Y.ravel()
    p = np.column_stack((x, y))
    v = da.values.ravel()
    idx = ~np.isnan(v)
    interpolator = LinearNDInterpolator(p[idx, :], v[idx])

    vv = interpolator(pv)

    return vv
