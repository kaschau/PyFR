from pyfr.multicomp.eos.base import BaseEOS
import itertools as it
import numpy as np


class cpgEOS(BaseEOS):
    name = 'cpg'
    def __init__(self, cfg):
        super().__init__(cfg)

        self.input_props = [
            'MW',
            'cp0',
        ]

    def compute_consts(self, props, consts):
        self.consts = consts
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

    def diff_con_to_pri(self, cons, diff_cons):
        consts = self.consts
        cp0 = consts["cp0"]
        MW = consts["MW"]
        Ru = consts["Ru"]
        ns = consts['ns']
        ndims = len(cons) - (ns - 1) - 2

        rhoYk = cons[0:ns]
        *rhouvw, rhoE = cons[ns::]
        diff_rhoY= diff_cons[0:ns]
        *diff_rhouvw, diff_rhoE = diff_cons[ns::]

        rho = sum(rhoYk)
        diff_rho = sum(diff_rhoY)

        # Compute primiatives
        pris = self.con_to_pri(cons)
        p = pris[0]
        uvw = pris[1:ndims+1]
        T = pris[ndims+1]

        # Divide rhoY by ρ
        Yk = [rhoY / rho for rhoY in rhoYk]

        # Compute mixture properties
        Rmix = 0.0
        cp = 0.0
        for n, Y in enumerate(Yk):
            Rmix += Y / MW[n]
            cp += Y * cp0[n]
        Rmix *= Ru

        # Compute the temperature, pressure
        e = rhoE / rho - 0.5 * sum(v * v for v in uvw)

        # Velocity gradients: ∂u⃗ = 1/ρ·[∂(ρu⃗) - u⃗·∂ρ]
        diff_uvw = [(diff_rhov - v*diff_rho) / rho
                    for diff_rhov, v in zip(diff_rhouvw, uvw)]

        # Species gradients: ∂Y⃗ = 1/ρ·[∂(ρY⃗) - Y⃗·∂ρ]
        diff_Yk = [(diff_rhoY - Y*diff_rho) / rho
                    for diff_rhoY, Y in zip(diff_rhoY, Yk)]

        # Build temperature gradient
        diff_T = 1.0/rho*(diff_rhoE - rhoE/rho*diff_rho) - sum([i*j for i,j in zip(uvw,diff_uvw)])

        for n, (Y, diff_Y) in enumerate(zip(Yk, diff_Yk)):
            e_Y =  T*(cp0[n] - Ru/MW[n])
            diff_T -= e_Y*diff_Y
        diff_T /= (cp - Rmix)

        # Build pressure gradient
        diff_p = Rmix*T*diff_rho + rho*Rmix*diff_T + rho*T*Ru*sum([dY/M for dY,M in zip(diff_Yk,MW)])

        return [diff_p, *diff_uvw, diff_T, *diff_Yk[0:-1]]