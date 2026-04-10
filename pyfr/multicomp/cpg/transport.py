from pyfr.multicomp.cpg.eos import CPGEos
from pyfr.multicomp.cpg.species import CPGTransportSpecies


class ConstantTransport(CPGEos):
    species_cls = CPGTransportSpecies

    def __init__(self, cfg, *args, **kwargs):
        super().__init__(cfg, *args, **kwargs)
        self.trans = 'constant-props'
