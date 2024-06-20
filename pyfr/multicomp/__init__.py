from pyfr.multicomp import eos, transport
import numpy as np
from pathlib import Path
import yaml

def complete_species(key, usersp, refsp):
    # A function to collect the data in order of species listed in the input spdata
    # returns a numpy array of the data.
    prop = []
    for sp in usersp.keys():
        try:
            prop.append(usersp[sp][key])
        except KeyError:
            try:
                prop.append(refsp[sp][key])
            except KeyError:
                raise KeyError(
                    f"You want to use species {sp}, but did not provide a {key}, and it is not in the PyFR species database."
                )
        except TypeError:
            try:
                prop.append(refsp[sp][key])
            except TypeError:
                raise TypeError(
                    "The top level in your spieces data input yaml file must only be species names."
                )
            except KeyError:
                raise KeyError(
                    f"You want to use species {sp}, but did not provide a {key}, and it is not in the PyFR species database."
                )
    if isinstance(prop[0], str):
        return prop
    else:
        return np.array(prop, dtype=np.float64)

def find_species_input(file_or_list):
    is_file = False

    test_paths = [
        # Assume we are given a relative path
        Path.cwd(),
        # Otherwise, search our database for the file
        Path(__file__).parent / 'database',
    ]
    for path in test_paths:
        fullpath = path / file_or_list
        if Path(fullpath).exists():
            is_file = True
            break

    if is_file:
        with open(fullpath, "r") as f:
            usersp = yaml.load(f, Loader=yaml.SafeLoader)
    else:
        usersp = {key: dict() for key in file_or_list.replace(' ','').split(",")}

    return usersp