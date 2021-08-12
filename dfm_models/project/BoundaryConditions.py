"""Object representing FM boundary conditions

cody.l.johnson@erdc.dren.mil

"""

import re

import pandas as pd

from dfm_models._internal import validate_file


class Boundary:

    """Docstring for Boundary."""

    def __init__(self):
        """TODO: to be defined.

        :Docstring for Boundary.: TODO

        """


class WaterLevelBC(Boundary):

    """Docstring for WaterLevelBC."""

    def __init__(self):
        """TODO: to be defined.

        :Docstring for WaterLevelBC.: TODO

        """


class AstroWaterLevelBC(WaterLevelBC):

    """Docstring for WaterLevelBC."""

    def __init__(self, pli_fn, bc_fn):
        """TODO: to be defined.

        :Docstring for WaterLevelBC.: TODO

        """

        self.pli_fn = validate_file(pli_fn)
        self.boundary_points = parse_pli(self.pli_fn)

        self.bc_fn = validate_file(bc_fn)
        self.get_harmonic_constituents()

    def get_harmonic_constituents(self):
        """TODO: Docstring for get_harmonic_constatns.

        :function: TODO
        :returns: TODO

        """
        harcons = {}

        for point_id in self.boundary_points.point_id:

            harcons[point_id] = parse_bc(self.bc_fn, point_id)

        self.harcons = harcons


#################
#   functions   #
#################
def parse_pli(pli_fn):
    """TODO: Docstring for parse_pli.

    :function: TODO
    :returns: TODO

    """
    boundary_points = pd.read_csv(
        pli_fn, skiprows=2, sep="\s+", names=["x", "y", "point_id"]  # noqa: W605
    )
    return boundary_points


def parse_bc(bc_fn, point_id, quantity="astronomic"):
    """Parse bc

    :function: TODO
    :returns: TODO

    """

    # regex to find record for point_id
    record = re.compile(point_id)

    # loop over lines in *.bc file to find record
    with open(bc_fn, mode="r") as f:

        i = 0
        for line in f:

            # conditionally execute
            if record.search(line):
                break

            # increment line index until a match on record
            i += 1

    # particular amount of header lines and data for each record
    skiprows = i + 8
    harcons = pd.read_csv(
        bc_fn,
        sep="\s+",  # noqa: W605
        names=["component", "amplitude", "phase"],
        skiprows=skiprows,
        nrows=37,
    )

    return harcons
