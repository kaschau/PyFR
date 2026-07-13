import numpy as np

from pyfr.multicomp.chem.reaction import Reaction
from pyfr.util import subclass_where


class Chemistry:
    def __init__(self, cfg, *args, **kwargs):
        super().__init__(cfg, *args, **kwargs)

        self.reactions = []
        for i, rxn_sect in enumerate(self._reaction_sects):
            rtype = rxn_sect['rtype']
            rxn_cls = subclass_where(Reaction, name=rtype)
            self.reactions.append(rxn_cls(i, rxn_sect))

        # Populate each reaction with species data for self-contained expressions
        nu_f = self.nu_f
        nu_b = self.nu_b
        aij = self.aij
        MWs = self.MWs
        for rxn in self.reactions:
            rxn.populate(self.ns, MWs, nu_f[rxn.index], nu_b[rxn.index],
                         aij[rxn.index], self.sp_idx)

    @property
    def nr(self):
        return len(self.reactions)

    def reactions_by_type(self, name):
        return [rxn for rxn in self.reactions if rxn.name == name]

    @property
    def nu_f(self):
        mat = np.zeros((self.nr, self.ns))
        for rxn in self.reactions:
            for name, coeff in rxn.reactants.items():
                mat[rxn.index, self.sp_idx(name)] = coeff
        return mat

    @property
    def nu_b(self):
        mat = np.zeros((self.nr, self.ns))
        for rxn in self.reactions:
            for name, coeff in rxn.products.items():
                mat[rxn.index, self.sp_idx(name)] = coeff
        return mat

    @property
    def aij(self):
        mat = np.ones((self.nr, self.ns))
        for rxn in self.reactions:
            if hasattr(rxn, 'efficiencies'):
                mat[rxn.index, :] = getattr(rxn, 'default_efficiency', 1.0)
                for name, eff in rxn.efficiencies.items():
                    mat[rxn.index, self.sp_idx(name)] = eff
        return mat
