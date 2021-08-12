"""Classes and functions for calibration.

cody.l.johnson@erdc.dren.mil

"""

import numpy as np


def write_correction_file(corr_fn, correction_factors):
    """write out correction factors to Delft3D-FM readable forcing file

    :function: TODO
    :returns: TODO

    """

    with open(corr_fn, "w") as f:

        for boundary_point in correction_factors.keys():
            f.write("[forcing]\n")
            f.write(f"Name                            = {boundary_point}\n")
            f.write("Function                        = astronomic-correction\n")
            f.write("Quantity                        = astronomic component\n")
            f.write("Unit                            = -\n")
            f.write("Quantity                        = waterlevelbnd amplitude\n")
            f.write("Quantity                        = waterlevelbnd phase\n")
            f.write("Unit                            = deg\n")

            for component, corr_fac in correction_factors[boundary_point].items():

                f.write(f"{component:<8s}{corr_fac:>8.6f} 0.0\n")

            f.write("\n")


def compute_phase_offset(BC, const_errors, stations, thresh=0.01):
    """compute amplitude corrections for AstroWaterLevelBC instance BC

    :function: TODO
    :returns: TODO

    """

    phase_offsets = {}

    # iterate over boundary_points in BC
    for boundary_point in BC.boundary_points.itertuples():

        point_id = boundary_point.point_id
        harcons = BC.harcons[point_id]
        harcons = harcons[harcons.amplitude > thresh]

        phase_offsets[point_id] = {}

        # loop over harmonic constituents
        for component in harcons.component:

            eps = idw(boundary_point, const_errors, stations, component, p=2)

            phase_offsets[point_id][component] = eps

    return phase_offsets


def compute_amp_correction_factor(BC, const_errors, stations, thresh=0.01):
    """compute amplitude corrections for AstroWaterLevelBC instance BC

    :function: TODO
    :returns: TODO

    """

    correction_factors = {}

    # iterate over boundary_points in BC
    for boundary_point in BC.boundary_points.itertuples():

        point_id = boundary_point.point_id
        harcons = BC.harcons[point_id]
        harcons = harcons[harcons.amplitude > thresh]

        correction_factors[point_id] = {}

        # loop over harmonic constituents
        for component in harcons.component:

            eps = idw(boundary_point, const_errors, stations, component)

            a = harcons[harcons.component == component]["amplitude"].values[0]

            # calculate correction factor
            a_p = a - eps
            A = a / a_p

            correction_factors[point_id][component] = A

    return correction_factors


def idw(boundary_point, const_errors, stations, component, p=1):
    """Calculate inverse distance weighted correction quanitty (amplitude or phase) for harmonic component

    @params:

    boundary_point = NamedTuple from Boundary object containing a boudnary_point attribute
    const_errors = Dictionary containing station_name to harmonic analysis error Series
    stations = Dictionary containing station_name to COOPs objects
    component = Harmonic constant of concern

    """

    errors = []
    weights = []

    bm_lon = boundary_point.x
    bm_lat = boundary_point.y

    bm_ll = (bm_lat, bm_lon)

    for station in stations.keys():

        sn_lat = stations[station].lat
        sn_lon = stations[station].lon
        sn_ll = (sn_lat, sn_lon)

        # make sure that components match
        try:
            e = const_errors[station].loc[component]
        except KeyError:
            raise KeyError(
                f"Component: {component} was not found in error table for {station}"
            )

        w = weight_function(bm_ll, sn_ll, p)

        errors.append(e)
        weights.append(w)

    errors = np.array(errors)
    weights = np.array(weights)

    correction = np.sum(weights * errors) / np.sum(weights)

    return correction


def weight_function(bm_ll, sn_ll, p):
    """calculate weight function for boundary point and observation point

    inputs:
        bm_ll = lat,lon tuple of boundary point
        sn_ll = lat,lon  ytuple of obs point
        p = exponent in weight function

    """

    lat1 = bm_ll[0]
    lon1 = bm_ll[1]

    lat2 = sn_ll[0]
    lon2 = sn_ll[1]

    d = haversine_distance(lat1, lon1, lat2, lon2)

    return 1 / (d ** p)


def haversine_distance(lat1, lon1, lat2, lon2):
    """calculate distance using Haversine formula

    Assumes xy are in degrees longitude, latitude in WGS

    """

    lat1 *= np.pi / 180
    lon1 *= np.pi / 180

    lat2 *= np.pi / 180
    lon2 *= np.pi / 180

    R = 6371e3  # radius of earth in meters

    DelLat = lat2 - lat1
    DelLon = lon2 - lon1

    a = np.sin(DelLat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(DelLon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return R * c
