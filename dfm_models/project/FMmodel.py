"""Object representing FM model

cody.l.johnson@erdc.dren.mil

"""

import re
from collections.abc import Iterable

import numpy as np
import pandas as pd
import pytide
import xarray as xr

from dfm_models._internal import validate_file, validate_project


class Model:

    """Top-level object for organizing FM model input/output"""

    def __init__(self, name, project_dir, results=None, **kwargs):
        """TODO: to be defined.

        :Docstring for model.: TODO

        """

        self.name = name
        self.project_dir = validate_project(project_dir)
        self.mdu = find_MDU(project_dir)
        self.params = load_MDU(self.mdu)

        # Check type of results
        if results is not None:
            if not isinstance(results, Results):
                raise TypeError(
                    "Results object of incorrect type. Initialize as FMmodel.Results class"
                )

        # save results if None or Results object
        self.results = results


class Grid:

    """Docstring for Grid."""

    def __init__(self):
        """TODO: to be defined.

        :Docstring for Grid.: TODO

        """


class UnstructGrid(Grid):

    """Docstring for UnstructGrid."""

    def __init__(self, fn):
        """TODO: to be defined.

        :Docstring for UnstructGrid.: TODO

        """
        self.fn = fn
        self.data = self.load_ugrid()

    def load_ugrid(self):
        data = xr.open_dataset(self.fn)
        return data


class Results:

    """Class to organize various results from an FM model"""

    def __init__(self, his_output):

        self.his_output = his_output


class hisOutput:

    """Class representation of history file.

    This is intented to package pre/post-processing code to work with
    history data.

    Add data processing methods here.

    """

    def __init__(self, his_fn):

        self.his_fn = validate_file(his_fn)
        self.data = load_his(self.his_fn)

    def water_level_comp_data(self, stations):
        """Create dictionary of water level time series and their demeaned
        values.

        :function: TODO
        :returns: TODO

        """

        water_level = {}
        water_level_LMSL = {}

        # make sure station is iterable
        if not isinstance(stations, Iterable):
            raise TypeError(
                "Station list is not iterable. Check that argument contains a list of valid station names."
            )

        # Check if any stations are in his_output
        if not any(item in stations for item in self.data.station_name):

            print(f"No stations in {stations} are in history output.")
            self.water_level = water_level
            self.water_level_LMSL = water_level_LMSL

        else:

            for station_name in stations:

                try:

                    data = (
                        self.data["waterlevel"]
                        .sel(station_name=station_name)
                        .to_dataframe()["waterlevel"]
                    )
                    data.rename()
                    water_level[station_name] = data

                    LMSL = data.mean()
                    water_level_LMSL[station_name] = data - LMSL

                except KeyError:
                    print(
                        f"{station_name} from station list doesn't match with any history output station names."
                    )

            self.water_level = water_level
            self.water_level_LMSL = water_level_LMSL

    def harmonic_analysis(
        self,
        consts=[
            "K1",
            "O1",
            "Sa",
            "Ssa",
            "P1",
            "Q1",
            "M2",
            "OO1",
            "M11",
            "S2",
            "J1",
            "Rho1",
            "K2",
            "S1",
            "Mm",
            "MSf",
            "N2",
            "MS4",
        ],
    ):
        """perform harmonic analysis on water level data, optionally using constants in consts.

        :function: TODO
        :returns: TODO

        """

        harcons = {}

        for station_name in self.data["station_name"].values:

            data = self.data["waterlevel"].sel(station_name=station_name).values
            time = self.data.time.values

            # use full set of constants if consts isn't defined
            if consts is None:
                wt = pytide.WaveTable()
                consts = wt.known_constituents()
            else:
                wt = pytide.WaveTable(consts)

            f, vu = wt.compute_nodal_modulations(time)
            w = wt.harmonic_analysis(data, f, vu)

            amplitude = np.abs(w)
            phase = np.angle(w, deg=True)

            harcons[station_name] = pd.DataFrame(
                {"amplitude": amplitude, "phase": phase},
                index=[const.upper() for const in consts],
            )

        self.harcons = harcons


#################
#   functions   #
#################
def find_MDU(project_dir):
    """glob for mdu in project_dir

    :function: TODO
    :returns: TODO

    """

    mdu = list(project_dir.glob("*.mdu"))[0]

    if not mdu.exists():
        print(f"No MDU found in {project_dir}.")
        print("Check MDU file existence.")

    return mdu


def find_net(project_dir):
    """glob for mdu in project_dir

    :function: TODO
    :returns: TODO

    """

    net_nc = list(project_dir.glob("*_net.nc"))[0]

    if not net_nc.exists():
        print(f"No Net file found in {project_dir}.")
        print("Check net file existence.")

    return net_nc


def load_MDU(mdu):
    """parse MDU and store parameter key/values in dictionary

    :function: TODO
    :returns: TODO

    """

    params = {}
    key_value = re.compile("(.*)=(.*)#.*")

    with open(mdu, "r") as f:
        lines = f.readlines()

        for line in lines:
            match = key_value.match(line)

            if match is None:
                continue

            else:

                key = match.group(1).strip()
                value = match.group(2).strip()
                params[key] = value

    return params


def load_his(his_fn):
    """load history file into Xarray DataSet w/ pre-processing

    :function: TODO
    :returns: TODO

    """

    _his = xr.open_dataset(his_fn)
    _his["station_name"] = _his["station_name"].str.decode("utf-8")

    return _his.swap_dims({"stations": "station_name"}).drop_vars(["station_id"])
