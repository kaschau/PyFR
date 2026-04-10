import numpy as np

from pyfr.multicomp import RU
from pyfr.multicomp.base import BaseEos
from pyfr.multicomp.tpg.fitting import (
    fit_poly, eval_nasa7_cp, eval_nasa7_h, eval_nasa7_s,
    eval_nasa9_cp, eval_nasa9_h, eval_nasa9_s,
    eval_over_ranges
)
from pyfr.multicomp.tpg.species import TPGSpecies


def _eval_cp_R(sp, T):
    """Evaluate cp/R for a species at temperature(s) T using active data."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ranges = sp.thermo_ranges
    coeffs = sp.thermo_coeffs

    if sp.eval_mode == 'fast':
        # Optimized: single range, variable-degree poly (last 2 are h,s consts)
        c = coeffs[0]
        nc = len(c) - 2
        result = np.zeros_like(T)
        for i in range(nc):
            result += c[i] * T**i
        return result
    else:
        # Strict: use raw NASA evaluation
        if sp._raw_model == 'nasa7':
            return eval_over_ranges(T, ranges, coeffs, eval_nasa7_cp)
        else:
            return eval_over_ranges(T, ranges, coeffs, eval_nasa9_cp)


def _eval_h_R(sp, T):
    """Evaluate h/R for a species at temperature(s) T using active data."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ranges = sp.thermo_ranges
    coeffs = sp.thermo_coeffs

    if sp.eval_mode == 'fast':
        c = coeffs[0]
        nc = len(c) - 2
        h_const = c[-2]
        result = np.zeros_like(T)
        for i in range(nc):
            result += c[i] * T**(i + 1) / (i + 1)
        return result + h_const
    else:
        if sp._raw_model == 'nasa7':
            return eval_over_ranges(T, ranges, coeffs, eval_nasa7_h)
        else:
            return eval_over_ranges(T, ranges, coeffs, eval_nasa9_h)


class TPGEos(BaseEos):
    name = 'tpg'
    species_cls = TPGSpecies

    def __init__(self, cfg, *args, **kwargs):
        super().__init__(cfg, *args, **kwargs)

        self.fit_tol = cfg.getfloat('multi-component', 'fit-tol', 0.0)
        self.T_min = (cfg.getfloat('multi-component', 'T-min')
                      if cfg.hasopt('multi-component', 'T-min') else None)
        self.T_max = (cfg.getfloat('multi-component', 'T-max')
                      if cfg.hasopt('multi-component', 'T-max') else None)
        if self.fit_tol > 0:
            self._fit_thermo()
            self.eval_mode = 'fast'
            self.T_iter_count = 4
            self.T_iter_method = 'halley'
            for sp in self.species:
                sp.bind_fast()
        else:
            self.eval_mode = 'strict'
            self.T_iter_count = 7
            self.T_iter_method = 'newton'
            for sp in self.species:
                sp.bind_strict()

    def _fit_thermo(self):
        for sp in self.species:
            fit_lo = self.T_min if self.T_min is not None else sp._raw_ranges[0]
            fit_hi = self.T_max if self.T_max is not None else sp._raw_ranges[-1]
            Ts = np.linspace(fit_lo, fit_hi, 500)

            cp_fn = eval_nasa7_cp if sp._raw_model == 'nasa7' else eval_nasa9_cp
            h_fn = eval_nasa7_h if sp._raw_model == 'nasa7' else eval_nasa9_h
            s_fn = eval_nasa7_s if sp._raw_model == 'nasa7' else eval_nasa9_s

            cp = eval_over_ranges(Ts, sp._raw_ranges, sp._raw_coeffs, cp_fn)
            h_ref = eval_over_ranges(Ts, sp._raw_ranges, sp._raw_coeffs, h_fn)
            s_ref = eval_over_ranges(Ts, sp._raw_ranges, sp._raw_coeffs, s_fn)

            cp_coeffs = fit_poly(Ts, cp, self.fit_tol)

            h_fit = sum(c * Ts**(i + 1) / (i + 1)
                        for i, c in enumerate(cp_coeffs))
            h_const = float(np.mean(h_ref - h_fit))

            s_fit = cp_coeffs[0] * np.log(Ts)
            for i in range(1, len(cp_coeffs)):
                s_fit += cp_coeffs[i] * Ts**i / i
            s_const = float(np.mean(s_ref - s_fit))

            sp.thermo_ranges = [fit_lo, fit_hi]
            sp.thermo_coeffs = [cp_coeffs + [h_const, s_const]]

    def validate_Y(self, Yk):
        for Y in Yk:
            if np.any(Y < 0) or np.any(Y > 1):
                raise ValueError('Species mass fraction out of range [0, 1]')

    def pri_to_con(self, pris):
        ns = self.ns
        ndims = len(pris) - (ns - 1) - 2

        p = np.asarray(pris[0])
        T = np.asarray(pris[ndims + 1])
        vs = pris[1:ndims + 1]
        Yk = list(pris[ndims + 2:])
        Yns = 1.0 - sum(Yk)
        Yk.append(Yns)

        self.validate_Y(Yk)

        # Mixture R and enthalpy
        Rmix = 0.0
        h = 0.0
        for n, Y in enumerate(Yk):
            sp = self.species[n]
            Rmix += Y / sp.MW
            h += Y * _eval_h_R(sp, T) * RU / sp.MW
        Rmix *= RU

        # Density
        rho = p / (Rmix * T)

        # Momentum
        rhovs = [rho * v for v in vs]

        # Total energy: rhoE = rho*h - p + 0.5*rho*|v|^2
        rhoE = rho * h - p + 0.5 * rho * sum(v * v for v in vs)

        # Species mass
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

        # Mixture R
        Rmix = sum(Y / sp.MW for Y, sp in zip(Yk, self.species)) * RU

        # Internal energy
        e = rhoE / rho - 0.5 * sum(v * v for v in vs)

        # Newton iteration for temperature
        T_lo = max(sp._raw_ranges[0] for sp in self.species)
        T_hi = min(sp._raw_ranges[-1] for sp in self.species)
        T = np.ones(np.asarray(e).shape) * 0.5 * (T_lo + T_hi)
        for _ in range(10):
            h = 0.0
            cp = 0.0
            for n, Y in enumerate(Yk):
                sp = self.species[n]
                cp += Y * _eval_cp_R(sp, T) * RU / sp.MW
                h += Y * _eval_h_R(sp, T) * RU / sp.MW
            error = e - (h - Rmix * T)
            T -= error / (-cp + Rmix)

        p = rho * Rmix * T

        return [p, *vs, T, *Yk[:-1]]

    def diff_con_to_pri(self, cons, diff_cons):
        ns = self.ns
        ndims = len(cons) - ns - 1

        rhoYk = cons[:ns]
        *rhouvw, rhoE = cons[ns:]
        diff_rhoY = diff_cons[:ns]
        *diff_rhouvw, diff_rhoE = diff_cons[ns:]

        rho = sum(rhoYk)
        diff_rho = sum(diff_rhoY)

        # Compute primitives
        pris = self.con_to_pri(cons)
        p = pris[0]
        uvw = pris[1:ndims + 1]
        T = pris[ndims + 1]
        Yk = [rhoY / rho for rhoY in rhoYk]

        e = rhoE / rho - 0.5 * sum(v * v for v in uvw)

        # Velocity gradients
        diff_uvw = [(diff_rhov - v * diff_rho) / rho
                    for diff_rhov, v in zip(diff_rhouvw, uvw)]

        # Species gradients
        diff_Yk = [(diff_rhoY - Y * diff_rho) / rho
                   for diff_rhoY, Y in zip(diff_rhoY, Yk)]

        # Temperature gradient
        diff_T = (1.0 / rho * (diff_rhoE - rhoE / rho * diff_rho)
                  - sum(u * du for u, du in zip(uvw, diff_uvw)))

        Rmix = 0.0
        cp = 0.0
        MW = np.array([sp.MW for sp in self.species])
        for n, (Y, diff_Y) in enumerate(zip(Yk, diff_Yk)):
            sp = self.species[n]
            Rmix += Y / sp.MW

            cp_sp = _eval_cp_R(sp, T) * RU / sp.MW
            cp += cp_sp * Y

            hk = _eval_h_R(sp, T) * RU / sp.MW
            e_Y = hk - T * RU / sp.MW
            diff_T -= e_Y * diff_Y

        Rmix *= RU
        diff_T /= (cp - Rmix)

        # Pressure gradient
        diff_p = (Rmix * T * diff_rho + rho * Rmix * diff_T
                  + rho * T * RU * sum(dY / M for dY, M
                                       in zip(diff_Yk, MW)))

        return [diff_p, *diff_uvw, diff_T, *diff_Yk[:-1]]
