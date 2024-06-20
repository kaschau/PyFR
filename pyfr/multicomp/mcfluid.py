from pyfr.util import subclass_where
from pyfr.multicomp import complete_species, find_species_input
from pyfr.multicomp.eos import BaseEOS
from pyfr.multicomp.transport import BaseTransport
from pathlib import Path
import yaml

class MCFluid:
    def __init__(self, cfg):

        self.eos = cfg.get('multi-component','eos')
        self.trans = cfg.get('multi-component','transport', None)

        eos_data = subclass_where(BaseEOS, name=self.eos)()
        if self.trans != 'None':
            trans_data = subclass_where(BaseTransport, name=self.trans)()

        # Save the prims <-> cons functions
        self.pri_to_con = eos_data.pri_to_con

        # Merge the lists of required data
        self.input_props = eos_data.input_props
        if self.trans != 'None':
            self.input_props |= trans_data.input_props

        # Merge the constants lists
        self.consts = eos_data.consts
        if self.trans != 'None':
            self.consts |= trans_data.consts
        self.consts['Ru'] = 8314.46261815324
        self.consts['avogadro'] = 6.02214076e+26
        self.consts['kb'] = 1.380649e-23
        self.consts['epsilon0'] = 8.854187812773345e-12

        # Get our species names
        file_or_list = cfg.get('multi-component', 'species')
        usersp = find_species_input(file_or_list)

        # HACK: Default to unity lewis
        if self.trans == 'constant-props':
            for key in usersp:
                if 'Le' not in usersp[key].keys():
                    usersp[key]['Le'] = 1.0
        # HACK: Default to unity lewis

        self.consts['ns'] = len(usersp)
        self.consts['names'] = [key for key in usersp]

        # Now load reference species
        relpath = str(Path(__file__).parent)
        with open(f"{relpath}/database/species_library.yaml", "r") as f:
            refsp = yaml.load(f, Loader=yaml.SafeLoader)

        # Now fill in all the property data
        for key in self.input_props.keys():
            self.input_props[key] = complete_species(key, usersp, refsp)

        # Now we can compute/fill in any constants
        eos_data.compute_consts(self.input_props, self.consts)
        if self.trans != 'None':
            trans_data.compute_consts(self.input_props, self.consts)

    @staticmethod
    def get_species_names(cfg):
        file_or_list = cfg.get('multi-component', 'species')
        return list(find_species_input(file_or_list).keys())
