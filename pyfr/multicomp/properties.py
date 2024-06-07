import numpy as np

from pathlib import Path
import yaml

def get_yaml(fname):

    # Assume we are given a relative path
    fullpath = Path.cwd() / fname
    if Path(fullpath).exists():
        return fullpath
    # Otherwise, search our database for the file
    dbasepath = Path(__file__).parent / 'database' / fname
    if Path(dbasepath).exists():
        return dbasepath
    # Finally, assume we were passed a list of species
    return False


class BaseProperties:
    def __init__(self, cfg):
        self.eos = cfg.get("multi-component","eos")
        self.data = dict()
        self.data["Ru"] = 8314.46261815324

    @staticmethod
    def get_species_names(cfg):

        file_or_list = cfg.get("multi-component","species")
        isfile = get_yaml(file_or_list)
        if isfile:
            pass
        else:
            species_names = file_or_list.split(",")

        return species_names

    @staticmethod
    def get_num_species(cfg):

        file_or_list = cfg.get("multi-component","species")
        isfile = get_yaml(file_or_list)
        if isfile:
            pass
        else:
            species_names = file_or_list.split(",")

        return len(species_names)

    def fill_data(self):
        #TODO How are we filling in all the data
        self.data['MWinv'][0] = 1.0/10.0
        self.data['MWinv'][1] = 1.0/11.0
        self.data['MWinv'][2] = 1.0/12.0
        self.data['MWinv'][3] = 1.0/13.0
        self.data['cp0'][0] = 1000.0
        self.data['cp0'][1] = 1001.0
        self.data['cp0'][2] = 1002.0
        self.data['cp0'][3] = 1003.0


class ThermoProperties(BaseProperties):
    name = 'thermo'
    def __init__(self, cfg):
        super().__init__(cfg)

        self.species_names = self.get_species_names(cfg)
        self.ns = len(self.species_names)
        self.data["MW"] = np.zeros(self.ns)
        self.data["MWinv"] = np.zeros(self.ns)
        if self.eos == "cpg":
            self.data["cp0"] = np.zeros(self.ns)
        elif self.eos == "tpg":
            self.data["nasa7"] = np.zeros((self.ns, 15))

        self.fill_data()

class TransportProperties(ThermoProperties):
    name = 'transport'
    def __init__(self, cfg):
        super().__init__()

        if self.transport == "constant":
            self.data["mu0"] = np.zeros(self.ns)
            self.data["kappa0"] = np.zeros(self.ns)
            self.data["Le"] = np.zeros(self.ns)
        elif self.eos == "kinetic-theory":
            self.data["muPoly"] = np.zeros(self.ns,5)
            self.data["kappaPoly"] = np.zeros(self.ns,5)
            self.data["dijPoly"] = np.zeros(self.ns,5)
