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
    def validate_data(consts):
        pass

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
        Ru = consts['Ru']
        MW = consts['MW']
        N7 = consts['NASA7'] * Ru/MW[:, np.newaxis]
        for n, Y in enumerate(pris[ndims+2::]+[Yns]):
            m = np.where(T <= N7[n,0], 8, 1)
            h += (  T*(N7[n, m + 0]
                   + T*(N7[n, m + 1] / 2.0
                   + T*(N7[n, m + 2] / 3.0
                   + T*(N7[n, m + 3] / 4.0
                   + T*(N7[n, m + 4] / 5.0))))) + N7[n, m + 5]) * Y

        # Compute density
        rho = p/(Rmix*T)

        # Multiply velocity components by rho
        rhovs = [rho * c for c in pris[1 : ndims + 1]]

        # Compute the total energy
        rhok = 0.5 * rho * sum(c * c for c in pris[1 : ndims + 1])
        rhoe = rho * h - p
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

        N7 = consts['NASA7'] * consts['Ru'] / consts['MW'][:, np.newaxis]
        N7[:,0] = consts['NASA7'][:,0]
        # Iterate on T, start at 300K
        T = np.ones(rho.shape)*300.0
        error = np.ones(rho.shape)
        niter = 0
        tol = 1e-8
        while np.max(np.abs(error)) > tol and niter < 100:
            h = 0.0
            cp = 0.0
            for n, Y in enumerate(Yk+[Yns]):
                m = np.where(T <= N7[n,0], 8, 1)
                cp += (     N7[n, m + 0]
                       + T*(N7[n, m + 1]
                       + T*(N7[n, m + 2]
                       + T*(N7[n, m + 3]
                       + T*(N7[n, m + 4] ))))) * Y

                h += (  T*(N7[n, m + 0]
                       + T*(N7[n, m + 1] /2.0
                       + T*(N7[n, m + 2] /3.0
                       + T*(N7[n, m + 3] /4.0
                       + T*(N7[n, m + 4] /5.0)))))
                       +    N7[n, m + 5]) * Y
            error = e - (h - Rmix * T)
            # Newtons Method
            T = T - error / (-cp - Rmix)
            niter += 1
            # print(niter, np.max(T), np.max(np.abs(error)))

        p = rho*Rmix*T

        return [p, *vs, T, *Yk]