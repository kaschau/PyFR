from os import path

class BaseProperties:
    def __init__(self, path_or_list):

        if path.exists(path_or_list):
            pass
        else:
            self.species_names = path_or_list.split(",")

        self.ns = len(self.species_names)
        self.data = dict()

class ThermoProperties(BaseProperties):
    def __init__(self, path_or_list):
        super().__init__(path_or_list)

class TransportProperties(BaseProperties):
    def __init__(self, path_or_list):
        super().__init__()

class ChemistyProperties(ThermoProperties):
    def __init__(self, path_or_list):
        super().__init__()