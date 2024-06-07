from os import path

class BaseProperties:
    def __init__(self, cfg):

        self.species_names = self.get_species_names(cfg)
        self.ns = len(self.species_names)
        self.data = dict()

    @staticmethod
    def get_species_names(cfg):

        path_or_list = cfg.get("multi-component","species")
        if path.exists(path_or_list):
            pass
        else:
            species_names = path_or_list.split(",")

        return species_names

    @staticmethod
    def get_num_species(cfg):

        path_or_list = cfg.get("multi-component","species")
        if path.exists(path_or_list):
            pass
        else:
            species_names = path_or_list.split(",")

        return len(species_names)



class ThermoProperties(BaseProperties):
    def __init__(self, path_or_list):
        super().__init__(path_or_list)

class TransportProperties(BaseProperties):
    def __init__(self, path_or_list):
        super().__init__()

class ChemistyProperties(ThermoProperties):
    def __init__(self, path_or_list):
        super().__init__()