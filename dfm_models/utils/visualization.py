"""Visualization tools for comparing model results to observations

cody.l.johnson@erdc.dren.mil

"""

import cartopy.crs as ccrs
import holoviews as hv
import matplotlib.pyplot as plt
import numpy as np


def spatial_stat(
    xs,
    ys,
    stat,
    labels,
    extent=[-93.5, -87, 28, 31],
    quantity_str="normalized rmse [% tidal range]",
    s=80,
    cmap="turbo",
    vmin=0,
    vmax=30,
    tbuff=-0.05,
):
    fig = plt.figure(figsize=(16, 6.5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.coastlines(resolution="10m", color="black", linewidth=1)
    im = ax.scatter(
        xs,
        ys,
        s=s,
        c=stat,
        transform=ccrs.PlateCarree(),
        zorder=10,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    for x, y, l in zip(xs, ys, labels):
        ax.text(
            x=x + tbuff,
            y=y + tbuff,
            s=l,
            horizontalalignment="right",
            transform=ccrs.PlateCarree(),
            zorder=101,
        )

    cbar = plt.colorbar(im)
    cbar.set_label(quantity_str)

    return fig, ax, cbar


def one2one(obs, mod, quantity_str="water level [m, NAVD88]", lims=None, ax=None):
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
