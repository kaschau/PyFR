from pyfr.multicomp.eos.base import BaseEOS
import itertools as it
from collections import deque
import numpy as np
from scipy.optimize import minimize


class tpgEOS(BaseEOS):
    name = 'tpg'
    def __init__(self, cfg):
        super().__init__(cfg)

        self.input_props = {
            'MW': None,
            'NASA7': None,
        }

        self.consts = {
            'MW': None,
            'NASA7': None,
        }

    def compute_consts(self):
        self.consts['MW'] = self.input_props['MW']

        # Determine if we are using strict property cals or fast
        prop_calc = self.cfg.get('multi-component', 'property-calc', 'strict')
        if prop_calc == 'strict':
            self.consts['NASA7'] = self.input_props['NASA7']
        else:
            N7 = self.input_props['NASA7']
            Tmin = self.cfg.getfloat('multi-component', 'T-min', 300.0)
            Tmax = self.cfg.getfloat('multi-component', 'T-max', 5000.0)
            Ts = np.linspace(Tmin, Tmax, 500)

            self.consts['NASA7'] = np.empty((self.consts['ns'], 7))
            for n in range(self.consts['ns']):
                m = np.where(Ts <= N7[n, 0], 8, 1)

                # Fit h
                h = (  Ts*(N7[n, m + 0]
                     + Ts*(N7[n, m + 1] / 2.0
                     + Ts*(N7[n, m + 2] / 3.0
                     + Ts*(N7[n, m + 3] / 4.0
                     + Ts*(N7[n, m + 4] / 5.0))))) + N7[n, m + 5])
                h_poly = np.polynomial.Polynomial.fit(Ts, h, 5)

                coeffs = h_poly.convert().coef
                coeffs *= [1, 1, 2, 3, 4, 5]
                # reorder so a0 is first
                dque = deque(coeffs)
                dque.rotate(-1)
                coeffs = list(dque)

                # solve for a6 (entorpy)
                s = (  np.log(Ts)*N7[n, m + 0]
                            +(Ts *(N7[n, m + 1]
                            + Ts *(N7[n, m + 2] / 2.0
                            + Ts *(N7[n, m + 3] / 3.0
                            + Ts *(N7[n, m + 4] / 4.0))))) + N7[n, m + 6])

                s_new = (  np.log(Ts)*coeffs[0]
                            +(Ts *(coeffs[1]
                            + Ts *(coeffs[2] / 2.0
                            + Ts *(coeffs[3] / 3.0
                            + Ts *(coeffs[4] / 4.0))))))

                def f(a6):
                    return np.sum((s - (s_new + a6))**2)

                a6 = minimize(f, N7[n, 7]).x[0]

                coeffs.append(a6)
                self.consts['NASA7'][n, :] = coeffs


    def pri_to_con(self, pris):
        consts = self.consts
        ns = consts['ns']
        ndims = len(pris) - (ns - 1) - 2
        p, T = pris[0], pris[ndims + 1]

        # Compute ns species
        Yns = 1.0 - sum(pris[ndims+2::])

        # Check mass fractions all 0<Y<1
        self.validate_Y_ics(it.chain(pris[ndims+2::],[Yns]))

        # Compute mixture properties
        Rmix = 0.0
        for n, Y in enumerate(it.chain(pris[ndims+2::],[Yns])):
            Rmix += Y/consts['MW'][n]
        Rmix *= consts['Ru']

        # Compute h
        h = 0.0
        Ru = consts['Ru']
        MW = consts['MW']
        N7 = consts['NASA7'] * Ru/MW[:, np.newaxis]
        N7[:, 0] = consts['NASA7'][:, 0]
        for n, Y in enumerate(it.chain(pris[ndims+2::],[Yns])):
            if len(N7[n]) == 15: # strict
                m = np.where(T <= N7[n,0], 8, 1)
            elif len(N7[n]) == 7:
                m = 0
            else:
                raise ValueError("NASA7 Issue.")
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
        rhoYk = [rho * c for c in it.chain(pris[ndims+2::],[Yns])]

        return [*rhoYk, *rhovs, rhoE]

    def con_to_pri(self, cons):
        consts = self.consts
        ns = consts['ns']
        ndims = len(cons)-(ns-1)-2

        rhoY = cons[0 : ns]
        rho = sum(rhoY)
        rhoE = cons[-1]
        # Divide momentum components by rho
        vs = [rhov / rho for rhov in cons[ns : ns + ndims]]

        # Species Mass Fraction
        Yk = [rhoYk / rho for rhoYk in rhoY]

        # Compute mixture properties
        Rmix = 0.0
        for n, Y in enumerate(Yk):
            Rmix += Y / consts["MW"][n]
        Rmix *= consts["Ru"]

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
            for n, Y in enumerate(Yk):
                if N7[n].shape[0] == 15: # strict
                    m = np.where(T <= N7[n,0], 8, 1)
                elif N7[n].shape == 7:
                    m = 0
                else:
                    raise ValueError("NASA7 Issue.")
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

        return [p, *vs, T, *Yk[0:-1]]