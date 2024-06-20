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
        p, T = pris[0], pris[ndims + 1]

        # Compute ns species
        Yns = 1.0 - sum(pris[ndims+2::])
        # Compute mixture properties
        Rmix = 0.0
        for n, Y in enumerate(pris[ndims+2::]+[Yns]):
            Rmix += Y/consts['MW'][n]
        Rmix *= consts['Ru']

        # Compute h
        h = 0.0
        Tinv = 1.0 / T
        To2 = T / 2.0
        T2o3 = T**2 / 3.0
        T3o4 = T**3 / 4.0
        T4o5 = T**4 / 5.0
        N7 = consts['NASA7']
        Ru = consts['Ru']
        MW = consts['MW']
        for n, Y in enumerate(pris[ndims+2::]+[Yns]):
            m = np.where(T <= N7[n,0], 8, 1)
            hns = (N7[n,m+0] + N7[n,m+1]*To2
                             + N7[n,m+2]*T2o3
                             + N7[n,m+3]*T3o4
                             + N7[n,m+4]*T4o5
                             + N7[n,m+5]*Tinv) * T * Ru/MW[n]
            h += hns * Y
        # Compute density
        rho = p/(Rmix*T)

        # Multiply velocity components by rho
        rhovs = [rho * c for c in pris[1 : ndims + 1]]

        # Compute the total energy
        rhok = 0.5 * rho * sum(c * c for c in pris[1 : ndims + 1])
        #HERE
        rhoe = rho * h - p
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
        Yns = 1.0 - sum(Yk)
        # Compute mixture properties
        Rmix = 0.0
        for k,Y in enumerate(Yk+[Yns]):
            Rmix += Y/consts['MW'][k]
        Rmix *= consts['Ru']

        # Internal energu
        e = rhoE/rho - 0.5 * sum(v * v for v in vs)

        N7 = consts['NASA7']
        Ru = consts['Ru']
        MW = consts['MW']
        # Iterate on T, start at 300K
        T = 300.0
        for _ in range(20):
            h = 0.0
            cp = 0.0
            Tinv = 1.0 / T
            To2 = T / 2.0
            T2 = T*T
            T3 = T2*T
            T4 = T3*T
            T2o3 = T2 / 3.0
            T3o4 = T3 / 4.0
            T4o5 = T4 / 5.0
            for n, Y in range(ns):
                m = 8 if T <= N7[n,0] else 1
                cps = (N7[n,m+0] + N7[n,m+1]*T
                                 + N7[n,m+2]*T2
                                 + N7[n,m+3]*T3
                                 + N7[n,m+4]*T4) * T * Ru/MW[n]
                hs = (N7[n,m+0] + N7[n,m+1]*To2
                                + N7[n,m+2]*T2o3
                                + N7[n,m+3]*T3o4
                                + N7[n,m+4]*T4o5
                                + N7[n,m+5]*Tinv) * T * Ru/MW[n]
                cp += cps * Y
                h += hs * Y
            error = e - (h - Rmix * T)
            T = T - error / (-cp - Rmix)

        assert np.all(np.abs(error) < 1e-8)

        p = rho*Rmix*T

        return [p, *vs, T, *Yk]