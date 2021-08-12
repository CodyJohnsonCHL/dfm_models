import pathlib
from datetime import datetime as dt


def validate_COOPs_loaded(Obs):
    """Simple validation that COOPs stations are properly loaded"""
    try:
        if not bool(Obs.COOPs):
            raise TypeError(
                "Load COOPs stations before downloading predicted water levels."
            )
    except AttributeError:
        raise AttributeError("Load COOPs stations before downloading data.")


def validate_harcon(container):
    """ensure harcons exists and are correct type

    :function: TODO
    :returns: TODO

    """

    try:
        if not bool(container.harcons):
            raise TypeError(
                "harcons must be a dictionary mapping station names to DataFrame of harmonic constiuents."
            )
    except AttributeError:
        raise AttributeError(
            f"{container} has no harcons. Check that harmonic analysis has been completed."
        )


def validate_project(project_dir):
    """validate path to FM model project dir

    :function: various check on top level DFM project directory
    :returns: PurePath object of DFM project

    """

    if not isinstance(project_dir, pathlib.Path):

        project_dir = pathlib.Path(project_dir)

        if not project_dir.is_dir():

            print(f"{project_dir} is not a valid directory.")
            print("Check path specification.")
            raise

    return project_dir


def validate_file(file_path):
    """check type and validate his_fn. Raise error is incorrect

    :function: TODO
    :returns: TODO

    """

    if not isinstance(file_path, pathlib.Path):

        file_path = pathlib.Path(file_path)

        if not file_path.is_file():

            raise TypeError(
                f"{file_path} is not a valid file.Check path specification."
            )

    return file_path


def validate_datetime(datetime):
    if not isinstance(datetime, dt):
        raise TypeError(f"{datetime} should be of type datetime.")


def validate_variable(variables):
    valid_vars = ["salinity", "water_temp", "surf_el", "water_u", "water_v"]

    if not set(variables) <= set(valid_vars):

        raise TypeError(f"{variables} should be a subset of {valid_vars}")
