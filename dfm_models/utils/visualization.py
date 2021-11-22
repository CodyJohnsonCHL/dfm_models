"""Visualization tools for comparing model results to observations

cody.l.johnson@erdc.dren.mil

"""

import geoviews as gv
import geoviews.feature as gf
import geoviews.tile_sources as gts
import holoviews as hv
import matplotlib.pyplot as plt
import numpy as np
from bokeh.models import HoverTool
from geoviews import opts


def spatial_stat(stats, fn):
    # data
    lon = hv.Dimension("lon", label="Longitude [deg]")
    lat = hv.Dimension("lat", label="Latitude [deg]")
    hv_stats = hv.Table(stats, kdims=[lon, lat])
    cols = stats.columns[2:-2].values
    clabels = {
        "nrmse": "normalized RMSE [% range]",
        "nrmse_tide": "normalized RMSE [% tidal range]",
        "rmse_cm": "RMSE [cm]",
        "r2": "r-squared [.]",
    }

    # hover tool
    tooltips = [
        ("Station", "@station"),
        ("# obs.", "@number"),
        ("Normalized RMSE [% range]", "@nrmse"),
        ("Normalized RMSE [% tidal range]", "@nrmse_tide"),
        ("RMSE [cm]", "@rmse_cm"),
        ("r-squared", "@r2"),
    ]
    hover = HoverTool(tooltips=tooltips)

    # style
    psize = 10
    cst_lw = 1.25

    # Holoviews options
    cOpts = opts.LineContours(line_width=cst_lw)
    overOpts = opts.Overlay(aspect=6.5 / 3, responsive=True)

    # generate holomap
    holomap = hv.HoloMap(kdims="Statistic")
    for col in cols:

        clabel = clabels[col]  # colorbar text label

        # options for points
        pOpts = opts.Points(
            size=psize,
            color=col,
            cmap="turbo",
            colorbar=True,
            clabel=clabel,
            tools=[hover],
            clim=(0, hv_stats[col].max()),
        )

        # put together
        overlay = (
            gf.coastline(scale="10m").opts(cOpts)
            * gv.Points(hv_stats).opts(pOpts)
            * gts.EsriImagery
        )

        # map
        holomap[col] = overlay.opts(overOpts)

    # save output
    gv.save(holomap, fn)

    return holomap


def one2one(obs, mod, quantity_str="water level [m, MSL]", lims=None, ax=None):
    std = obs.std()
    std_shift = np.sin(np.pi / 4) * std

    if lims is None:
        lims = [np.min([obs.min(), mod.min()]), np.max([obs.max(), mod.max()])]

    lower_bound_x = [lims[0] + std_shift, lims[1]]
    lower_bound_y = [lims[0], lims[1] - std_shift]
    upper_bound_x = [lims[0], lims[1] - std_shift]
    upper_bound_y = [lims[0] + std_shift, lims[1]]

    fill_between_x = lims
    fill_between_y1 = [lims[0] - std_shift, lims[1] - std_shift]
    fill_between_y2 = [lims[0] + std_shift, lims[1] + std_shift]

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(obs, mod, 15)
        ax.plot(lims, lims, color="k", linestyle="--")

        ax.plot(lower_bound_x, lower_bound_y, color="k", lw=0.75)
        ax.plot(upper_bound_x, upper_bound_y, color="k", lw=0.75)
        ax.fill_between(
            fill_between_x, fill_between_y1, fill_between_y2, alpha=0.2, color="gray"
        )
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel(f"observed {quantity_str}")
        ax.set_ylabel(f"modeled {quantity_str}")

        return fig, ax
    else:
        ax.scatter(obs, mod, 15)
        ax.plot(lims, lims, color="k", linestyle="--")

        ax.plot(lower_bound_x, lower_bound_y, color="k", lw=0.75)
        ax.plot(upper_bound_x, upper_bound_y, color="k", lw=0.75)
        ax.fill_between(
            fill_between_x, fill_between_y1, fill_between_y2, alpha=0.2, color="gray"
        )
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel(f"observed {quantity_str}")
        ax.set_ylabel(f"modeled {quantity_str}")
        return ax


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
