import numpy as np

from pyfr.multicomp import KB, RU, EPS0
from pyfr.multicomp.tpg.eos import TPGEos
from pyfr.multicomp.tpg.fitting import (
    fit_poly, eval_nasa7_cp, eval_nasa9_cp,
    eval_over_ranges, omega22_at, astar_at
)
from pyfr.multicomp.tpg.species import TPGTransportSpecies


class KineticTheory(TPGEos):
    species_cls = TPGTransportSpecies

    def __init__(self, cfg, *args, **kwargs):
        super().__init__(cfg, *args, **kwargs)
        self.trans = 'kinetic-theory'
        self._fit_transport()

    def _fit_transport(self):
        tol = self.fit_tol

        # Global temperature range for transport fitting
        # Default: intersection of all species' thermo ranges
        T_lo = (self.T_min if self.T_min is not None
                else max(sp._raw_ranges[0] for sp in self.species))
        T_hi = (self.T_max if self.T_max is not None
                else min(sp._raw_ranges[-1] for sp in self.species))

        for sp in self.species:
            Ts = np.linspace(T_lo, T_hi, 50)
            logTs = np.log(Ts)
            sqrtTs = np.sqrt(Ts)

            # Point-wise collision integral evaluation (matches Cantera)
            omga22 = np.empty_like(Ts)
            omga11 = np.empty_like(Ts)
            for i, T in enumerate(Ts):
                ts = T * KB / sp.well_depth
                o22 = omega22_at(ts, sp.deltastar_self)
                ast = astar_at(ts, sp.deltastar_self)
                omga22[i] = o22
                omga11[i] = o22 / ast

            mu = ((5.0/16.0) * np.sqrt(np.pi * sp.mass * KB * Ts)
                  / (np.pi * sp.diameter**2 * omga22))

            r_mass = sp.mass / 2.0
            diff_self = ((3.0/16.0) * np.sqrt(2.0*np.pi / r_mass)
                         * (KB*Ts)**1.5
                         / (np.pi * sp.diameter**2 * omga11))

            f_int = sp.MW / (RU * Ts) * diff_self / mu
            Tstar = Ts * KB / sp.well_depth
            Tstar_298 = KB * 298.0 / sp.well_depth
            fz_298 = (1.0 + np.pi**1.5 / np.sqrt(Tstar_298)
                      * (0.5 + 1.0/Tstar_298)
                      + (0.25*np.pi**2 + 2) / Tstar_298)
            fz_T = (1.0 + np.pi**1.5 / np.sqrt(Tstar)
                    * (0.5 + 1.0/Tstar)
                    + (0.25*np.pi**2 + 2) / Tstar)
            A_factor = 2.5 - f_int
            B_factor = (sp.rot_relax * fz_298/fz_T
                        + 2.0/np.pi * (5.0/3.0*sp.rot_dof + f_int))
            c1 = 2.0/np.pi * A_factor / B_factor

            cp_fn = (eval_nasa7_cp if sp._raw_model == 'nasa7'
                     else eval_nasa9_cp)
            cp_R = eval_over_ranges(Ts, sp._raw_ranges, sp._raw_coeffs,
                                    cp_fn)
            cv_int = cp_R - 2.5 - sp.rot_dof
            f_trans = 2.5 * (1.0 - c1*sp.rot_dof/1.5)
            f_rot = f_int * (1.0 + c1)
            kappa = (mu / sp.MW * RU
                     * (f_trans*1.5 + f_rot*sp.rot_dof + f_int*cv_int))

            mu_scaled = np.sqrt(mu / sqrtTs)
            sp.muPoly = fit_poly(logTs, mu_scaled, tol,
                                 w=1.0/np.abs(mu_scaled))

            kappa_scaled = kappa / sqrtTs
            sp.kappaPoly = fit_poly(logTs, kappa_scaled, tol,
                                    w=1.0/np.abs(kappa_scaled))

        self._fit_Dij(tol, T_lo, T_hi)

    def _fit_Dij(self, tol, T_lo, T_hi):
        ns = self.ns
        species = self.species

        Ts = np.linspace(T_lo, T_hi, 50)
        logTs = np.log(Ts)

        polys = [[None]*ns for _ in range(ns)]

        for n in range(ns):
            for n2 in range(n, ns):
                sp1, sp2 = species[n], species[n2]

                r_mass = sp1.mass * sp2.mass / (sp1.mass + sp2.mass)
                r_diam = 0.5 * (sp1.diameter + sp2.diameter)
                r_well = np.sqrt(sp1.well_depth * sp2.well_depth)
                r_dipole = np.sqrt(sp1.dipole * sp2.dipole)

                r_deltastar = (
                    0.5 * r_dipole**2
                    / (4*np.pi * EPS0 * r_well * r_diam**3)
                )

                if sp1.polar == sp2.polar:
                    f_well, f_diam = 1.0, 1.0
                else:
                    kp = sp1 if sp1.polar else sp2
                    knp = sp2 if sp1.polar else sp1
                    d3np = knp.diameter**3
                    d3p = kp.diameter**3
                    alphastar = knp.polarizability / d3np
                    dipolestar = kp.dipole / np.sqrt(
                        4*np.pi * EPS0 * d3p * kp.well_depth
                    )
                    xi = (1.0 + 0.25 * alphastar * dipolestar**2
                          * np.sqrt(kp.well_depth / knp.well_depth))
                    f_well = xi**2
                    f_diam = xi**(-1/6)

                r_well *= f_well
                r_diam *= f_diam

                diff = np.empty_like(Ts)
                for i, T in enumerate(Ts):
                    Tstar = T * KB / r_well
                    o22 = omega22_at(Tstar, r_deltastar)
                    ast = astar_at(Tstar, r_deltastar)
                    o11 = o22 / ast

                    diff[i] = ((3.0/16.0)
                               * np.sqrt(2.0*np.pi / r_mass)
                               * (KB * T)**1.5
                               / (np.pi * r_diam**2 * o11))

                inv_diff_scaled = 1.0 / (diff / Ts**1.5)
                w = 1.0 / np.abs(inv_diff_scaled)
                poly = fit_poly(logTs, inv_diff_scaled, tol, w=w)

                polys[n][n2] = poly
                polys[n2][n] = poly

        for n in range(ns):
            species[n].DijPoly = polys[n]
