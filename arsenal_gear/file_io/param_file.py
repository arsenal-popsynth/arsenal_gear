"""
param_file
==========

Functions for loading and checking the arsenal param file.
"""

import yaml

param_types = {
    "outputs": str,
    "formation": {
        "binaries": bool,
        "ages": {
            "type": str,
            "tform": float,
        },
        "metals": {
            "type": str,
            "metals": float,
        },
        "imf": {
            "type": str,
            "min_mass": float,
            "max_mass": float,
            "N": int,
        },
    },
    "evolution": {
        "SN": {
            "type": str,
            "min_mass": float,
            "max_mass": float,
            "calculate": [str],
            "energy": float,
        },
    },
}


def load_param_file(filename: str) -> dict:
    """
    Load a YAML parameter file and return its contents as a nested set of
    dictionaries.

    Parameters
    ----------
    filename : str
        The path to the YAML file.

    Returns
    -------
    dict
        The contents of the YAML file as a dictionary.
    """
    with open(filename, encoding="utf-8") as file:
        params = yaml.safe_load(file)
    return params


def check_param_file(
    params: dict,
) -> bool:  # pylint: disable=too-many-return-statements
    """
    Check if the parameter file is valid.

    Parameters
    ----------
    params : dict
        The parameters to check.

    Returns
    -------
    bool
        True if the parameter file is valid, False otherwise.
    """
    for k in params:
        if k not in param_types:
            print(f"Unknown parameter: {k}")
            return False
        if isinstance(param_types[k], dict):
            for sub_k in params[k]:
                if sub_k not in param_types[k]:
                    print(f"Unknown parameter: {k}.{sub_k}")
                    return False
                if not isinstance(params[k][sub_k], param_types[k][sub_k]):
                    print(
                        f"Incorrect type for {k}.{sub_k}: expected "
                        f"{param_types[k][sub_k]}, got {type(params[k][sub_k])}"
                    )
                    return False
        elif isinstance(param_types[k], list):
            if not isinstance(params[k], list):
                print(f"Incorrect type for {k}: expected list, got {type(params[k])}")
                return False
            for item in params[k]:
                if not isinstance(item, param_types[k][0]):
                    print(
                        f"Incorrect type for item in {k}: expected "
                        f"{param_types[k][0]}, got {type(item)}"
                    )
                    return False
        elif not isinstance(params[k], param_types[k]):
            print(
                f"Incorrect type for {k}: expected {param_types[k]}, got {type(params[k])}"
            )
            return False
    return True
