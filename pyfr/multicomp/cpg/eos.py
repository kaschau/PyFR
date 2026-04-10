import numpy as np

from pyfr.multicomp import RU
from pyfr.multicomp.base import BaseEos
from pyfr.multicomp.cpg.species import CPGSpecies


class CPGEos(BaseEos):
    name = 'cpg'
    species_cls = CPGSpecies

    def __init__(self, cfg, *args, **kwargs):
        super().__init__(cfg, *args, **kwargs)

    def pri_to_con(self, pris):
        ndims = len(pris) - (self.ns - 1) - 2
        p, T = np.asarray(pris[0]), np.asarray(pris[ndims + 1])
        vs = pris[1:ndims + 1]
        Yk = list(pris[ndims + 2:])
        Yns = 1.0 - sum(Yk)
        Yk.append(Yns)

        Rmix = sum(Y / sp.MW for Y, sp in zip(Yk, self.species)) * RU
        cpmix = sum(Y * sp.cp0 for Y, sp in zip(Yk, self.species))

        rho = p / (Rmix * T)
        rhovs = [rho * v for v in vs]
        rhoE = rho * cpmix * T - p + 0.5 * rho * sum(v*v for v in vs)
        rhoYk = [rho * Y for Y in Yk]

        return [*rhoYk, *rhovs, rhoE]

    def con_to_pri(self, cons):
        ns = self.ns
        ndims = len(cons) - ns - 1

        rhoY = cons[:ns]
        rho = sum(rhoY)
        rhoE = cons[-1]
        vs = [rhov / rho for rhov in cons[ns:ns + ndims]]
        Yk = [rhoYk / rho for rhoYk in rhoY]

        Rmix = sum(Y / sp.MW for Y, sp in zip(Yk, self.species)) * RU
        cpmix = sum(Y * sp.cp0 for Y, sp in zip(Yk, self.species))

        e = rhoE / rho - 0.5 * sum(v*v for v in vs)
        T = e / (cpmix - Rmix)
        p = rho * Rmix * T

        return [p, *vs, T, *Yk[:-1]]
