from pyfr.multicomp.eos.baseEOS import BaseEOS
import numpy as np


class tpgEOS(BaseEOS):
    name = 'tpg'
    def __init__(self):
        super().__init__()

        self.input_props = {
            'MW': None,
            'NASA7': None,
        }

        self.consts = {
            'MW': None,
            'NASA7': None,
        }

    @staticmethod
    def compute_consts(props, consts):
        consts['MW'] = props['MW']
        consts['NASA7'] = props['NASA7']

    def pri_to_con(self, pris):
        consts = self.consts
        ns = consts['ns']
        ndims = len(pris) - (ns - 1) - 2

        # Compute ns species
        Yns = 1.0 - np.sum(pris[ndims+2::])
        # Compute mixture properties
        Rmix = 0.0
        cp = 0.0
        for k,Y in enumerate(pris[ndims+2::]+[Yns]):
            Rmix += Y/consts['MW'][k]
            cp += Y*consts['cp0'][k]
        Rmix *= consts['Ru']

        # Compute density
        p, T = pris[0], pris[ndims + 1]
        rho = p/(Rmix*T)

        # Multiply velocity components by rho
        rhovs = [rho * c for c in pris[1 : ndims + 1]]

        # Compute the total energy
        rhok = 0.5 * rho * np.sum(c * c for c in pris[1 : ndims + 1])
        #HERE
        rhoe = rho*T*(cp-Rmix)
        #HERE
        rhoE = rhoe + rhok

        # Species mass
        rhoYk = [rho * c for c in pris[ndims + 2::]]

        return [rho, *rhovs, rhoE, *rhoYk]

    def con_to_pri(self, cons):
        consts = self.consts
        ns = consts['ns']
        ndims = len(cons)-(ns-1)-2

        rho, rhoE = cons[0], cons[ndims + 1]

        # Divide momentum components by rho
        vs = [rhov / rho for rhov in cons[1 : ndims + 1]]

        # Species Mass Fraction
        Yk = [rhoY / rho for rhoY in cons[ndims + 2 ::]]

        # Compute ns species
        Yns = 1.0 - np.sum(Yk)
        # Compute mixture properties
        Rmix = 0.0
        cp = 0.0
        for k,Y in enumerate(Yk+[Yns]):
            Rmix += Y/consts['MW'][k]
            cp += Y*consts['cp0'][k]
        Rmix *= consts['Ru']

        # Compute the temperature, pressure
        e = rhoE/rho - 0.5 * np.sum(v * v for v in vs)

        # Here
        T = e/(cp-Rmix)
        # Here

        p = rho*Rmix*T

        return [p, *vs, T, *Yk]