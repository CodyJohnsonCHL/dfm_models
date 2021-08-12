#!/usr/bin/env python
# coding: utf-8

# version 0.1   2020/01/08  -- Cody Johnson

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


# functions #
def write_bc_for_pli(bc_fn, gcm, pli, quantity, method, depth_avg):
    """
    append or write 3D boundary conditions for quantities

    bc_fn = path to write or append boundary condition data
    gcm = general circulation model output which contains boundary points (xr.DataArray)
    pli = table of boundary support points (pd.DataFrame)
    quantity = variable to output to the BC files (salinity or water_temp)
    depth_avg = flag to enable depth averaging
    """

    with open(bc_fn, "a") as f:

        gcm_refd, ref_date = assign_seconds_from_refdate(gcm)

        for _, (x_pli, y_pli, pli_point_name) in pli.iterrows():

            x_pli_east = x_pli + 360

            if (quantity == "salinity") or (quantity == "water_temp"):
                bc_data = interpolate_to_pli_point(
                    gcm_refd, quantity, x_pli_east, y_pli, pli_point_name, method
                )

                # check valid depth point
                if bc_data is None:
                    continue

                if depth_avg:
                    write_ts_record(f, bc_data, pli_point_name, quantity, ref_date)
                else:
                    write_t3d_record(f, bc_data, pli_point_name, quantity, ref_date)

            # in the case of velocity both components are interpolated
            elif quantity == "velocity":
                bc_data = interpolate_to_pli_point(
                    gcm_refd,
                    ["water_u", "water_v"],
                    x_pli_east,
                    y_pli,
                    pli_point_name,
                    method,
                )

                # check valid depth point
                if bc_data is None:
                    continue

                write_vector_t3d_record(f, bc_data, pli_point_name, quantity, ref_date)

            # in case of surf_el
            elif quantity == "surf_el":
                bc_data = interpolate_to_pli_point(
                    gcm_refd, quantity, x_pli_east, y_pli, pli_point_name, method
                )

                # check valid depth point
                if bc_data is None:
                    continue

                write_ts_record(f, bc_data, pli_point_name, quantity, ref_date)


def write_ts_record(f, bc_data, pli_point_name, quantity, ref_date):
    """
    append or write time series boundary conditions for depth averaged quantities

    f = file descriptor for writing boundary condition output
    bc_data = data at with percent bed coords (xr.DataFrame)
    pli_point_name = name for entry
    quantity = variable to output to the BC files
    ref_date = date used to calculate offset of the record in seconds
    """

    # get units for quantity
    if quantity == "salinity":
        quantbnd = "salinitybnd"
        units = "ppt"
    elif quantity == "water_temp":
        quantbnd = "temperaturebnd"
        units = "°C"
    elif quantity == "surf_el":
        quantbnd = "waterlevelbnd"
        units = "m"
    else:
        print('quantity needs to be either "salinity", "water_temp" or "surf_el"\n')
        raise ValueError

    # write a record
    f.write("[forcing]\n")
    f.write(f"Name                            = {pli_point_name}\n")
    f.write(f"Function                        = t3d\n")
    f.write(f"Time-interpolation              = linear\n")
    f.write(f"Vertical position type          = single\n")
    f.write(f"Vertical interpolation          = linear\n")
    f.write(f"Quantity                        = time\n")
    f.write(f"Unit                            = seconds since {ref_date}\n")
    f.write(f"Quantity                        = {quantbnd}\n")
    f.write(f"Unit                            = {units}\n")
    f.write(f"Vertical position               = 1\n")

    if quantity == "surf_el":
        for td, value in bc_data.to_dataframe()[quantity].iteritems():

            value = f"{value:.05f}"
            f.write(f"{td} {value}\n")

    else:
        # write data after converting to dataframe and iterating over the rows
        for td, values in bc_data.to_dataframe()[quantity].unstack(level=-1).iterrows():

            # take mean of values to get depth averaged
            value = values.mean()

            # see results of interpolation
            if value > 100.0:
                print(
                    f"Problem with {quantity} exceeding maximum allowed value: {values.max():.03f} ppt."
                )
            elif value < 0.0:
                print(
                    f"Problem with {quantity} becoming negative: {values.max():.03f} ppt."
                )
                print(f"Negative value for {quantity} has been set to 0.01 {units}.")
                value = 0.01

            value = f"{value:.05f}"
            f.write(f"{td} {value}\n")

    f.write("\n")


def write_vector_t3d_record(f, bc_data, pli_point_name, quantity, ref_date):
    """
    append or write 3D boundary conditions for quantities

    f = file descriptor for writing boundary condition output
    bc_data = data at with percent bed coords (xr.DataFrame)
    pli_point_name = name for entry
    quantity = variable to output to the BC files
    ref_date = date used to calculate offset of the record in seconds
    """

    if quantity == "velocity":
        vector = "uxuyadvectionvelocitybnd:ux,uy"
        quantbndx = "ux"
        quantbndy = "uy"
        x_comp = "water_u"
        y_comp = "water_v"
        units = "-"  # no units for velocity in example provided by Kees
    else:
        print('quantity should be "velocity"\n')
        raise ValueError

    # convert percent from bed into formated string
    pos_spec = [f"{perc:.02f}" for perc in bc_data.perc_from_bed.data]
    pos_spec_str = " ".join(pos_spec[::-1])  # reverse order for D3D

    # write a record
    f.write("[forcing]\n")
    f.write(f"Name                            = {pli_point_name}\n")
    f.write(f"Function                        = t3d\n")
    f.write(f"Time-interpolation              = linear\n")
    f.write(f"Vertical position type          = percentage from bed\n")
    f.write(f"Vertical position specification = {pos_spec_str}\n")
    f.write(f"Vertical interpolation          = linear\n")
    f.write(f"Quantity                        = time\n")
    f.write(f"Unit                            = seconds since {ref_date}\n")
    f.write(f"Vector                          = {vector}\n")

    # loop over number of vertical positions
    for vert_pos in range(1, len(pos_spec) + 1):
        f.write(f"Quantity                        = {quantbndx}\n")
        f.write(f"Unit                            = {units}\n")
        f.write(f"Vertical position               = {vert_pos}\n")

        f.write(f"Quantity                        = {quantbndy}\n")
        f.write(f"Unit                            = {units}\n")
        f.write(f"Vertical position               = {vert_pos}\n")

    # write data after converting to dataframe and iterating over the rows
    for td, values in (
        bc_data.to_dataframe()[[x_comp, y_comp]].unstack(level=0).iterrows()
    ):

        # get componets as array in order to format for d3d input
        x_comp_vals = values[x_comp].values[::-1]  # reverse order for D3D
        y_comp_vals = values[y_comp].values[::-1]  # reverse order for D3D
        values = [
            f"{x_comp_val:.03f} {y_comp_val:.03f}"
            for x_comp_val, y_comp_val in zip(x_comp_vals, y_comp_vals)
        ]
        values_str = " ".join(values)
        f.write(f"{td} {values_str}\n")

    f.write("\n")


def write_t3d_record(f, bc_data, pli_point_name, quantity, ref_date):
    """
    append or write 3D boundary conditions for quantities

    f = file descriptor for writing boundary condition output
    bc_data = data at with percent bed coords (xr.DataFrame)
    pli_point_name = name for entry
    quantity = variable to output to the BC files
    ref_date = date used to calculate offset of the record in seconds
    """

    # get units for quantity
    if quantity == "salinity":
        quantbnd = "salinitybnd"
        units = "ppt"
    elif quantity == "water_temp":
        quantbnd = "temperaturebnd"
        units = "°C"
    else:
        print('quantity needs to be either "salinity" or "water_temp"\n')
        raise ValueError

    # convert percent from bed into formated string
    pos_spec = [f"{perc:.02f}" for perc in bc_data.perc_from_bed.data]
    pos_spec_str = " ".join(pos_spec[::-1])  # reverse order for D3D

    # write a record
    f.write("[forcing]\n")
    f.write(f"Name                            = {pli_point_name}\n")
    f.write(f"Function                        = t3d\n")
    f.write(f"Time-interpolation              = linear\n")
    f.write(f"Vertical position type          = percentage from bed\n")
    f.write(f"Vertical position specification = {pos_spec_str}\n")
    f.write(f"Vertical interpolation          = linear\n")
    f.write(f"Quantity                        = time\n")
    f.write(f"Unit                            = seconds since {ref_date}\n")

    # loop over number of vertical positions
    for vert_pos in range(1, len(pos_spec) + 1):
        f.write(f"Quantity                        = {quantbnd}\n")
        f.write(f"Unit                            = {units}\n")
        f.write(f"Vertical position               = {vert_pos}\n")

    # write data after converting to dataframe and iterating over the rows
    for td, values in bc_data.to_dataframe()[quantity].unstack(level=-1).iterrows():

        # see results of interpolation
        if values.max() > 100.0:
            print(
                f"problem with {quantity} exceeding maximum allowed value: {values.max():.03f} ppt"
            )
        elif values.min() < 0.0:
            print(f"problem with {quantity} becoming negative: {values.max():.03f} ppt")
            print(f"Negative values for {quantity} has been set to 0.01 {units}.")
            values.where(values > 0.01, 0.01, inplace=True)

        values = [f"{value:.05f}" for value in values]
        values_str = " ".join(values[::-1])  # reverse order for D3D
        f.write(f"{td} {values_str}\n")

    f.write("\n")


def assign_seconds_from_refdate(gcm):
    """
    This func assigns seconds from a user specified ref date as coords.
    This is how D3D interpolates the boundary conditions in time.

    gcm = model output to add coords to
    """
    ref_date = gcm.time.data[0]
    ref_dt = pd.to_datetime(ref_date)
    ref_date_str = ref_dt.strftime("%Y-%m-%d %H:%M:%S")
    timedeltas = pd.to_datetime(gcm.time.data) - ref_dt
    seconds = timedeltas.days * 24 * 60 * 60 + timedeltas.seconds
    gcm = gcm.assign_coords(coords={"seconds_from_ref": ("time", seconds)})
    return gcm.swap_dims({"time": "seconds_from_ref"}), ref_date_str


def interpolate_to_pli_point(
    gcm_refd, quantity, x_pli_east, y_pli, pli_point_name, method
):
    """interpolates the quanitites to the sigma depths and pli coords

    gcm_refd = gcm with new time coordinates
    quantity = variable to output to the BC files (salinity or water_temp)
    x_pli_east = longitude of pli point in degrees east from meridian (GCM convention)
    y_pli = latitude
    """

    if quantity == "surf_el":

        # interpolate to pli point and drop data below bed level at nearest gcm_refd point
        bc_data = gcm_refd[quantity].interp(lon=x_pli_east, lat=y_pli, method=method)
        return bc_data

    else:

        # interpolate to pli point and drop data below bed level at nearest gcm_refd point
        bc_data = (
            gcm_refd[quantity]
            .interp(lon=x_pli_east, lat=y_pli, method=method)
            .dropna(dim="depth")
            .squeeze()
        )

        # add coordinate for percent from bed. D3D uses this in its bc file format
        try:
            gcm_refd_zb = bc_data.depth[-1]  # get bed level of gcm_refd point
        except IndexError:
            print(
                f"Depth invalid for {pli_point_name} at: {x_pli_east}, {y_pli}. Omitting point..."
            )
            return None

        perc_from_bed = 100 * (-1 * bc_data.depth + gcm_refd_zb) / gcm_refd_zb
        bc_data = bc_data.assign_coords(
            coords={"perc_from_bed": ("depth", perc_from_bed)}
        )

        return bc_data


### main loop ###
if __name__ == "__main__":

    ### arguments ###
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "nc",
        help="GCM NetCDF output containing boundary support points and duration of Delft3D simulation",
    )
    parser.add_argument(
        "quantity",
        help='GCM variable. Must be either "salintiy", "water_temp", or "velocity"',
    )
    parser.add_argument(
        "--pli-list",
        nargs="*",
        type=str,
        help="list of boundary support point polyline filenames",
        required=True,
        dest="pli_list",
    )
    parser.add_argument(
        "--bc-filename",
        help="Optional filename for Delft3D boundary condition filename",
        type=str,
        dest="bc_filename",
    )
    parser.add_argument(
        "--depth-avg",
        help="flag to enable depth averaged output",
        default=False,
        action="store_true",
        dest="depth_avg",
    )
    parser.add_argument(
        "--interp-method",
        help="flag to enable depth averaged output",
        default="linear",
        type=str,
        dest="method",
    )
    args = parser.parse_args()

    gcm = args.nc
    quantity = args.quantity
    pli_list = args.pli_list
    depth_avg = args.depth_avg
    method = args.method

    # validate arguments
    if (
        (quantity != "salinity")
        and (quantity != "water_temp")
        and (quantity != "velocity")
        and (quantity != "surf_el")
    ):
        print(
            f'<quantity> was specfied as {quantity}, but should be either "salinity" or "water_temp".'
        )
        raise ValueError

    # open gcm NetCDF output as Xarray dataset
    try:
        gcm = xr.open_dataset(Path(gcm), drop_variables="tau")
    except FileNotFoundError as e:
        print("<GCM output> should be path to GCM NetCDF output")
        raise e

    # set defualt boundary condition filename depending on quanitity
    bc_fn = args.bc_filename
    if bc_fn is None:
        if quantity == "salinity":
            bc_fn = Path("Salinity.bc")
        elif quantity == "water_temp":
            bc_fn = Path("Temperature.bc")
        elif quantity == "velocity":
            bc_fn = Path("Velocity.bc")
        elif quantity == "surf_el":
            bc_fn = Path("WaterLevel.bc")

    # pli files opened as Pandas DataFrames
    pli_points = []
    for pli_fn in pli_list:

        print(f"Reading in file: {pli_fn}")
        pli = pd.read_csv(
            pli_fn, sep="\s+", skiprows=2, header=None, names=["x", "y", "point_id"]
        )

        write_bc_for_pli(bc_fn, gcm, pli, quantity, method=method, depth_avg=depth_avg)

        # add points to list for visualization
        pli_points.append(pli)

    # concat pli points
    pli_points = pd.concat(pli_points)

    ### visualization ###
    # color map depending on quantity
    if quantity == "salinity":
        cmap = cmo.haline
    elif quantity == "water_temp":
        cmap = cmo.thermal
    elif quantity == "velocity":
        cmap = "jet"
    else:
        cmap = "jet"

    # setup orthographic projection for geographic data
    fig, ax = plt.subplots(
        1, 1, subplot_kw={"projection": ccrs.Orthographic(-91, 29)}, figsize=(16, 9)
    )

    # plot initial quantity at surface
    if (quantity == "salinity") or (quantity == "water_temp"):
        gcm[quantity].isel(time=0, depth=0).plot(
            ax=ax, transform=ccrs.PlateCarree(), cmap=cmap
        )
    elif quantity == "velocity":
        tmp = gcm.isel(time=0, depth=0)
        tmp["magnitude"] = np.sqrt(tmp["water_u"] ** 2 + tmp["water_v"] ** 2)
        tmp["magnitude"].plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            cbar_kwargs={"label": "velocity magnitude [m/s]"},
        )

    # add coastline for reference
    ax.add_feature(cfeature.COASTLINE, edgecolor="0.3")

    # boundary condition support points
    pli_points.plot.scatter(
        "x", "y", marker="x", color="k", ax=ax, transform=ccrs.PlateCarree()
    )

    fig.savefig("point_output_locations.png", bbox_inches="tight")
