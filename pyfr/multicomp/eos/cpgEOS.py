from pyfr.multicomp.eos.baseEOS import BaseEOS
import itertools as it
import numpy as np


class cpgEOS(BaseEOS):
    name = 'cpg'
    def __init__(self):
        super().__init__()

        self.input_props = {
            'MW': None,
            'cp0': None,
        }

        self.consts = {
            'MW': None,
            'cp0': None,
        }

    @staticmethod
    def validate_data(consts):
        pass

    @staticmethod
    def compute_consts(props, consts):
        consts['MW'] = props['MW']
        consts['cp0'] = props['cp0']

    def pri_to_con(self, pris):
        consts = self.consts
        ns = consts['ns']
        ndims = len(pris) - (ns - 1) - 2


        # Compute ns species
        Yns = 1.0 - sum(pris[ndims+2::])

        # Check mass fractions all 0<Y<1
        self.validate_Y_ics(it.chain(pris[ndims+2::],[Yns]))

        # Compute mixture properties
        Rmix = 0.0
        cp = 0.0
        for n,Y in enumerate(it.chain(pris[ndims+2::],[Yns])):
            Rmix += Y/consts['MW'][n]
            cp += Y*consts['cp0'][n]
        Rmix *= consts['Ru']

        # Compute density
        p, T = pris[0], pris[ndims + 1]
        rho = p/(Rmix*T)

        # Multiply velocity components by rho
        rhovs = [rho * c for c in pris[1 : ndims + 1]]

        # Compute the total energy
        rhok = 0.5 * rho * sum(c * c for c in pris[1 : ndims + 1])
        rhoe = rho*T*(cp-Rmix)
        rhoE = rhoe + rhok

        # Species mass
        rhoYk = [rho * c for c in it.chain(pris[ndims+2::],[Yns])]

        return [*rhoYk, *rhovs, rhoE]

    def con_to_pri(self, cons):
        consts = self.consts
        ns = consts['ns']
        ndims = len(cons) - (ns - 1) - 2

        rhoY = cons[0 : ns]
        rho = sum(rhoY)
        rhoE = cons[-1]

        # Divide momentum components by rho
        vs = [rhov / rho for rhov in cons[ns : ns + ndims]]

        # Species Mass Fraction
        Yk = [rhoYk / rho for rhoYk in rhoY]

        # Compute mixture properties
        Rmix = 0.0
        cp = 0.0
        for n, Y in enumerate(Yk):
            Rmix += Y / consts["MW"][n]
            cp += Y * consts["cp0"][n]
        Rmix *= consts["Ru"]

        # Compute the temperature, pressure
        e = rhoE / rho - 0.5 * sum(v * v for v in vs)
        T = e / (cp - Rmix)
        p = rho * Rmix * T

        return [p, *vs, T, *Yk[0:-1]]
