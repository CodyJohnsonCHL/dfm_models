"""Analysis tools for comparing model results to observations

cody.l.johnson@erdc.dren.mil

"""

import pandas as pd
import pytide
from numpy import abs, angle

from dfm_models._internal import validate_harcon


def compute_harcon_error(Observations, Results):
    """Calculate statistics between COOPs and FM model harmonic constiuents"""

    amplitude_error = {}
    phase_error = {}

    validate_harcon(Observations)
    validate_harcon(Results.his_output)

    common_stations = get_common_stations(
        Observations.harcons, Results.his_output.harcons
    )

    for station in common_stations:

        amplitude_error[station] = (
            (Results.his_output.harcons[station] - Observations.harcons[station])[
                "amplitude"
            ]
            .dropna()
            .sort_values(ascending=False)
        )

        phase_results = Results.his_output.harcons[station]["phase"]
        phase_obs = Observations.harcons[station]["phase"]
        common_comps = get_common_components(phase_results, phase_obs)

        errors = []

        for comp in common_comps:

            res = phase_results.loc[comp]
            obs = phase_obs.loc[comp]

            e = compute_phase_error(res, obs)
            errors.append(e)

        errors = pd.Series(errors, index=common_comps)
        errors.name = "phase"
        phase_error[station] = errors

    return amplitude_error, phase_error


def compute_phase_error(res, obs):
    """compute error in phase between 0 and 360"""

    if res < 0:
        res += 360

    if obs < 0:
        obs += 360

    if res >= obs:
        mu = res - obs
        return mu

    else:
        mu = res - obs
        mu += 360
        return mu


##########################
#    tidal harmonic      #
##########################
def harmonic_analysis(
    waterlevel,
    time,
    consts=[
        "K1",
        "O1",
        "P1",
        "M2",
        "Q1",
        "S2",
        "S1",
        "Mf",
        "N2",
        "K2",
        "J1",
        "Mm",
        "M4",
        "Sa",
        "Ssa",
    ],
):
    wt = pytide.WaveTable(consts)
    h = waterlevel.values
    f, vu = wt.compute_nodal_modulations(time)
    w = wt.harmonic_analysis(h, f, vu)
    hp = wt.tide_from_tide_series(time, w)
    return w, (h, hp, time), consts


def get_modulus_angle(w):
    modulus = abs(w)
    ang = angle(w, deg=True)
    return modulus, ang


##########################
#        getters        #
##########################
def get_common_stations(obs_harcons, res_harcons):
    """Get intersection of stations for the observations and results

    :function: TODO
    :returns: TODO

    """

    obs_stations = obs_harcons.keys()
    res_stations = res_harcons.keys()

    common_stations = list(set(obs_stations) & set(res_stations))

    if len(common_stations) == 0:
        print("There are no common stations between observations and results.")
        return None

    else:

        return common_stations


def get_common_components(results, obs):
    """get intersectin of constituents"""

    # cast to set for intersection operator
    res_comps = set(results.index)
    obs_comps = set(obs.index)

    common_comps = list(res_comps & obs_comps)
    return common_comps
