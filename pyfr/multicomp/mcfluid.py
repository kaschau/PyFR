from functools import cache

import numpy as np

from pyfr.multicomp import RU
from pyfr.multicomp.base import BaseEos
from pyfr.multicomp.chem import Chemistry
from pyfr.multicomp.cpg.transport import ConstantTransport
from pyfr.multicomp.tpg.transport import KineticTheory
from pyfr.multicomp.readers.cantera import read_cantera_yaml
from pyfr.util import subclass_where


class MCFluidBase:
    Ru = RU

    def __init__(self, cfg, *args, **kwargs):
        self.cfg = cfg
        self.eos = cfg.get('multi-component', 'eos')
        self.mixing_rule = cfg.get('multi-component', 'mixing-rule', 'Wilke')
        self.trans = 'none'

        filepath = cfg.get('multi-component', 'species')
        species_sects, reaction_sects = read_cantera_yaml(filepath)

        sp_cls = self.species_cls
        self.species = [
            sp_cls(i, name, spc_sect)
            for i, (name, spc_sect) in enumerate(species_sects.items())
        ]
        self._name_to_idx = {sp.name: sp.index for sp in self.species}
        self._reaction_sects = reaction_sects

    @staticmethod
    def get_species_names(cfg):
        filepath = cfg.get('multi-component', 'species')
        species_sects, _ = read_cantera_yaml(filepath)
        return list(species_sects.keys())

    @property
    def ns(self):
        return len(self.species)

    @property
    def sp_names(self):
        return [sp.name for sp in self.species]

    @property
    def MWs(self):
        return np.array([sp.MW for sp in self.species])

    def mcix(self, ndims):
        vix = self.ns
        Eix = self.ns + ndims
        rhoix = self.ns + ndims
        pix = rhoix + 1
        Tix = pix + 1
        return vix, Eix, rhoix, pix, Tix

    def sp_idx(self, name):
        return self._name_to_idx[name]

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.species[self._name_to_idx[key]]
        return self.species[key]

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def __reduce__(self):
        return (int, (id(self),))


@cache
def get_mcfluid(cfg, needs_transport=False):
    eos_name = cfg.get('multi-component', 'eos')
    chemistry = cfg.getbool('multi-component', 'chemistry', False)

    eos_cls = subclass_where(BaseEos, name=eos_name)
    if needs_transport:
        eos_cls = eos_cls.__subclasses__()[0]

    bases = []
    if chemistry:
        bases.append(Chemistry)
    bases.append(eos_cls)
    bases.append(MCFluidBase)

    cls = type('MCFluid', tuple(bases), {})
    return cls(cfg)
