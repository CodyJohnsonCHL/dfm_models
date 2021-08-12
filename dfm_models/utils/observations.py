"""Classes and functions for working with observational data.

cody.l.johnson@erdc.dren.mil

"""

import json
from collections.abc import Iterable
from tempfile import NamedTemporaryFile as tmp
from urllib.error import HTTPError
from urllib.request import urlretrieve

import pandas as pd

from dfm_models._internal import validate_COOPs_loaded


class Observations:

    """Docstring for Observations."""

    def __init__(self):
        """TODO: to be defined.

        :Docstring for Observations.: TODO

        """

    def load_COOPs_stations(self, stations):
        """TODO: Docstring for generate_water_level_comp.

        :function: TODO
        :returns: TODO

        """

        self.COOPs = {}

        # make sure station is iterable
        if not isinstance(stations, Iterable):
            raise TypeError(
                "Station list is not iterable. Check that argument contains a list of valid station names."
            )

        # make sure stations is a dict
        if not type(stations) is dict:
            raise TypeError(
                "Station should be dictionary mapping COOPs station names to station_ids."
            )

        for station_name, station_id in stations.items():

            # verify that COOPs station is valid
            request = f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{station_id}.json"
            try:
                _, _ = urlretrieve(request)
                station = COOPs(station_name, station_id)
                self.COOPs[station_name] = station
            except HTTPError:
                print(
                    f"{station_name} with id: {station_id} is not a valid COOPs station. Check request:"
                )
                print(request)
                pass

    def gen_COOPS_predicted_water_levels(self, begin_date, end_date, datum="MSL"):
        """TODO: Docstring for generate_water_level_comp.

        :function: TODO
        :returns: TODO

        """

        validate_COOPs_loaded(self)

        # differentiate object based on datum
        if datum == "MSL":

            self.predicted_water_levels_MSL = {}

            for station_name, station in self.COOPs.items():

                try:
                    station.download_prediction(datum, begin_date, end_date)
                    self.predicted_water_levels_MSL[
                        station_name
                    ] = station.water_level_prediction_MSL
                except HTTPError:
                    pass

        elif datum == "NAVD":

            self.predicted_water_level_NAVD = {}

            for station_name, station in self.COOPs.items():

                try:
                    station.download_prediction(datum, begin_date, end_date)
                    self.predicted_water_levels_NAVD[
                        station_name
                    ] = station.water_level_prediction_NAVD
                except HTTPError:
                    pass

        elif datum == "LMSL":

            self.predicted_water_levels_LMSL = {}

            for station_name, station in self.COOPs.items():

                try:
                    station.download_prediction("MSL", begin_date, end_date)
                    self.predicted_water_levels_LMSL[station_name] = (
                        station.water_level_prediction_MSL
                        - station.water_level_prediction_MSL.mean()
                    )
                except HTTPError:
                    pass

    def gen_COOPs_harcons(self):
        """TODO: Docstring for gen_COOPs_harcons.

        :function: TODO
        :returns: TODO

        """

        validate_COOPs_loaded(self)

        self.harcons = {}

        for station_name, station in self.COOPs.items():

            station.download_harcon()

            self.harcons[station_name] = station.harcon


class Station:

    """Docstring for Station."""

    def __init__(self):
        """TODO: to be defined.

        :Docstring for Station.: TODO

        """


class COOPs(Station):

    """Docstring for COOPs_Station."""

    def __init__(self, station_name, station_id):
        """TODO: to be defined.

        :Docstring for COOPs_Station.: TODO

        """
        self.station_name = station_name
        self.station_id = station_id

        self.get_station_metadata()

    def get_station_metadata(self):
        """TODO: Docstring for get_station_metadata.

        :function: TODO
        :returns: TODO

        """
        request = f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{self.station_id}.json"

        fn = tmp()

        try:
            txt, http = urlretrieve(request, fn.name)
        except HTTPError:
            print(
                f"Station metadata for {self.station_name}"
                f"with {self.station_id} was not found in CO-OPs database."
            )
            print(f"Check url for errors: {request}")
            raise

        metadata = json.load(fn)
        self.lat = metadata["stations"][0]["lat"]
        self.lon = metadata["stations"][0]["lng"]

    def download_prediction(self, datum, begin_date, end_date):
        """download predicted water level and append to attribute

        :function: TODO
        :returns: TODO

        """
        request = (
            f"https://tidesandcurrents.noaa.gov/api/datagetter?"
            f"begin_date={begin_date}"
            f"&end_date={end_date}"
            f"&station={self.station_id}"
            f"&product=predictions"
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
            print(
                f"{self.station_name} with {self.station_id} was not found in CO-OPs database."
            )
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

        data.name = "predicted_water_level"

        # differentiate object based on datum
        if datum == "MSL":
            self.water_level_prediction_MSL = data
        elif datum == "NAVD":
            self.water_level_prediction_NAVD = data

    def download_harcon(self):
        """TODO: Docstring for download_harcon.

        :function: TODO
        :returns: TODO

        """

        request = (
            f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/"
            f"{self.station_id}/harcon.json?units=metric"
        )

        fn = tmp()

        try:
            txt, http = urlretrieve(request, fn.name)
        except HTTPError:
            print(
                f"Harmonic consts for {self.station_name} with {self.station_id} was not found in CO-OPs database."
            )
            print(f"Check url for errors: {request}")
            raise

        harcon = json.load(fn)
        self.harcon = get_amp_phase(harcon)


##########################
#       functions        #
##########################
# parse json of harcons
def get_amp_phase(harcon):
    const = harcon["HarmonicConstituents"]
    names = []
    amps = []
    phases = []
    speeds = []

    for component in const:
        names.append(component["name"])
        amps.append(component["amplitude"])
        phases.append(component["phase_GMT"])
        speeds.append(component["speed"])
        pass

    return pd.DataFrame(
        index=names, data={"amplitude": amps, "phase": phases, "speed": speeds}
    )
