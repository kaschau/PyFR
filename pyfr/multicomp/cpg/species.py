from pyfr.multicomp.base import BaseSpecies


class CPGSpecies(BaseSpecies):
    def __init__(self, index, name, spc_sect):
        super().__init__(index, name, spc_sect)
        self.cp0 = float(spc_sect['cp0'])

    def cp_expr(self, var='T'):
        return str(self.cp0)

    def gbs_expr(self, var='T', var_logT='logT', var_Tinv='Tinv'):
        raise NotImplementedError('CPG does not support Gibbs computation')


class CPGTransportSpecies(CPGSpecies):
    def __init__(self, index, name, spc_sect):
        super().__init__(index, name, spc_sect)
        self.mu0 = float(spc_sect['mu0'])
        self.kappa0 = float(spc_sect['kappa0'])
        self.Le = float(spc_sect.get('Le', 1.0))
