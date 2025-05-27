from pyfr.util import subclass_where
from pyfr.multicomp import complete_species, find_species_input
from pyfr.multicomp.eos import BaseEOS
from pyfr.multicomp.transport import BaseTransport
from pathlib import Path
import yaml
import numpy as np

class MCFluid:
    def __init__(self, cfg, justTherm=False):

        self.cfg = cfg
        self.eos = cfg.get('multi-component','eos')
        if justTherm:
            self.trans = 'None'
        else:
            self.trans = cfg.get('multi-component','transport', 'None')

        eos_data = subclass_where(BaseEOS, name=self.eos)(cfg)
        if self.trans != 'None':
            trans_data = subclass_where(BaseTransport, name=self.trans)(cfg)

        # Save the prims <-> cons functions
        self.pri_to_con = eos_data.pri_to_con
        self.con_to_pri = eos_data.con_to_pri
        self.diff_con_to_pri = eos_data.diff_con_to_pri

        # Merge the lists of required data
        self.input_props = {k:None for k in eos_data.input_props}
        if self.trans != 'None':
            self.input_props |= {k:None for k in trans_data.input_props}

        # Get our species names
        file_or_list = cfg.get('multi-component', 'species')
        userdata = find_species_input(file_or_list)
        usersp = userdata['properties']

        # HACK: Default to unity lewis
        if self.trans == 'constant-props':
            for key in usersp:
                if 'Le' not in usersp[key].keys():
                    usersp[key]['Le'] = 1.0
        # HACK: Default to unity lewis

        # Now load reference species
        relpath = str(Path(__file__).parent)
        with open(f"{relpath}/database/species_library.yaml", "r") as f:
            refsp = yaml.load(f, Loader=yaml.SafeLoader)

        # Now fill in all the property data
        for key in self.input_props.keys():
            self.input_props[key] = complete_species(key, usersp, refsp)

        # Now we can compute/fill in constants
        self.consts = {}
        self.consts['Ru'] = 8314.46261815324
        self.consts['avogadro'] = 6.02214076e+26
        self.consts['kb'] = 1.380649e-23
        self.consts['epsilon0'] = 8.854187812773345e-12

        self.consts['ns'] = len(usersp)
        self.consts['names'] = [key for key in usersp]

        eos_data.compute_consts(self.input_props, self.consts)
        if self.trans != 'None':
            trans_data.compute_consts(self.input_props, self.consts)

        # Finally, merge reactions data to the consts, make them numpy arrays
        chem = cfg.getbool('multi-component','chemistry', False)
        if chem:
            dmap = {'single': np.float32, 'double':np.float64}
            dtype = dmap[cfg.get('backend','precision')]
            finfo = np.finfo(dtype)
            for k,v in userdata['reactions'].items():
                if not isinstance(v[0], str):
                    # clip by type
                    temp = np.array(v).clip(finfo.min, finfo.max)
                    self.consts[k] = temp
                else:
                    self.consts[k] = v

    @staticmethod
    def get_species_names(cfg):
        file_or_list = cfg.get('multi-component', 'species')
        return list(find_species_input(file_or_list)['properties'].keys())
