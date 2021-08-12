"""Visualization tools for comparing model results to observations

cody.l.johnson@erdc.dren.mil

"""

import holoviews as hv


def water_level_holomap(water_levels, label):
    """create a simple HoloMap with a key dim on station

    :function: TODO
    :returns: TODO

    """
    kdims = hv.Dimension("station")

    curves = {
        station: hv.Curve(wl, group=label, label=label)
        for station, wl in water_levels.items()
    }

    return hv.HoloMap(curves, kdims=kdims)


def harcon_error_holomap(harcon_errors, label, drop_small=True):
    """TODO: Docstring for harcon_error_holomap.

    :function: TODO
    :returns: TODO

    """

    kdims = hv.Dimension("station")
    bars = {}

    for station, error in harcon_errors.items():

        if drop_small:
            error = error[error >= 0.005]

        bars[station] = hv.Bars(error, group=label, label=label)

    return hv.HoloMap(bars, kdims=kdims)
