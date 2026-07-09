import math
import sys

import numpy as np

from pyfr.fluids.base import BaseFluid
from pyfr.fluids.constants import RU
from pyfr.fluids.nasa import (eval_nasa7_cp, eval_nasa7_h, eval_nasa7_s,
                              eval_over_ranges)
from pyfr.fluids.readers.cantera import read_cantera_yaml
from pyfr.fluids.species import CPGSpecies, TPGSpecies
from pyfr.fluids.transport import BaseTransport
from pyfr.util import subclass_where


class MCBaseFluid(BaseFluid):
    species_cls = None

    def _read_cfg(self, cfg):
        filepath = cfg.getpath('multi-component', 'species')
        species_sects, reaction_sects = read_cantera_yaml(filepath)

        # Constant properties may be given (or overridden) in the ini
        for name, sect in species_sects.items():
            for key in ('cp0', 'mu0', 'kappa0', 'Le'):
                v = cfg.get('multi-component', f'{key}-{name}', '')
                if v:
                    sect[key] = float(v)

        self.species = [
            self.species_cls(i, name, sect)
            for i, (name, sect) in enumerate(species_sects.items())
        ]
        self._species_sects = species_sects

        self.chemistry = cfg.getbool('multi-component', 'chemistry', False)
        if self.chemistry:
            self._read_reactions(reaction_sects)

    def _read_reactions(self, reaction_sects):
        from pyfr.fluids.reaction import Reaction

        if not all(hasattr(sp, 'gbs_expr') for sp in self.species):
            raise ValueError('chemistry requires NASA polynomial thermo; '
                             'set eos = mc-tpg')
        if not reaction_sects:
            raise ValueError('chemistry = true but the mechanism has no '
                             'reactions')

        names = self.sp_names
        ns, nr = self.ns, len(reaction_sects)

        self.reactions = []
        for i, rxn_sect in enumerate(reaction_sects):
            rxn_cls = subclass_where(Reaction, name=rxn_sect['rtype'])
            self.reactions.append(rxn_cls(i, rxn_sect))

        nu_f = np.zeros((nr, ns))
        nu_b = np.zeros((nr, ns))
        aij = np.ones((nr, ns))
        for rxn in self.reactions:
            for name, coeff in rxn.reactants.items():
                nu_f[rxn.index, names.index(name)] = coeff
            for name, coeff in rxn.products.items():
                nu_b[rxn.index, names.index(name)] = coeff
            for name, eff in getattr(rxn, 'efficiencies', {}).items():
                aij[rxn.index, names.index(name)] = eff

        MWs = [sp.MW for sp in self.species]
        for rxn in self.reactions:
            rxn.populate(ns, MWs, nu_f[rxn.index], nu_b[rxn.index],
                         aij[rxn.index], names.index)

        # Which species participate in any reaction
        self.species_participates = [
            bool(np.any(nu_f[:, n] != nu_b[:, n])) for n in range(ns)
        ]

    def reactions_by_type(self, name):
        return [rxn for rxn in self.reactions if rxn.name == name]

    @property
    def ns(self):
        return len(self.species)

    @property
    def sp_names(self):
        return [sp.name for sp in self.species]

    @property
    def nvars(self):
        return self.ns + self.ndims + 1

    @property
    def privars(self):
        vs = ['u', 'v', 'w'][:self.ndims]

        return ['p', *vs, 'T', *self.sp_names[:-1]]

    @property
    def convars(self):
        vs = ['rhou', 'rhov', 'rhow'][:self.ndims]

        return [*(f'rho{n}' for n in self.sp_names), *vs, 'E']

    @property
    def dualcoeffs(self):
        return self.convars

    @property
    def visvars(self):
        vs = ['u', 'v', 'w'][:self.ndims]
        varmap = {
            'pressure': ['p'],
            'velocity': vs,
            'temperature': ['T']
        }
        for sn in self.sp_names[:-1]:
            varmap[sn] = [sn]

        return varmap

    def _register_state(self):
        ndims, ns, nvars = self.ndims, self.ns, self.nvars
        MWs = [sp.MW for sp in self.species]
        reg = self._register_quantity

        reg('rho', (),
            lambda u: ' + '.join(f'{u}[{n}]' for n in range(ns)))
        reg('invrho', ('rho',), lambda u: '1.0/rho')
        reg('Y', ('invrho',),
            lambda u: [f'{u}[{n}]*invrho' for n in range(ns)])
        reg('v', ('invrho',),
            lambda u: [f'{u}[{ns + i}]*invrho' for i in range(ndims)])
        reg('E', (), lambda u: f'{u}[{nvars - 1}]')
        reg('R', ('Y',),
            lambda u: ' + '.join(f'Y[{n}]*{RU / MWs[n]}' for n in range(ns)))
        reg('e', ('E', 'rho', 'v', 'invrho'),
            lambda u: '(E - 0.5*rho*('
            + ' + '.join(f'(v[{i}])*(v[{i}])' for i in range(ndims))
            + '))*invrho')

        self._register_thermo()

    def pri_seed(self, pris):
        ndims = self.ndims
        Yk = list(pris[ndims + 2:])
        Yk.append(1.0 - sum(Yk))

        return {'p': pris[0], 'v': list(pris[1:ndims + 1]),
                'T': pris[ndims + 1], 'Y': Yk}


class MCCPGFluid(MCBaseFluid):
    name = 'mc-cpg'
    species_cls = CPGSpecies

    # Per-species enthalpy at an arbitrary temperature symbol
    def h_of_T(self, n, T):
        return f'{self.species[n].cp0}*({T})'

    def _register_thermo(self):
        ns = self.ns
        cp0s = [sp.cp0 for sp in self.species]
        reg = self._register_quantity

        reg('cpmix', ('Y',),
            lambda u: ' + '.join(f'Y[{n}]*{cp0s[n]}' for n in range(ns)))
        reg('T', ('e', 'cpmix', 'R'), lambda u: 'e/(cpmix - R)')
        reg('p', ('rho', 'R', 'T'), lambda u: 'rho*R*T')
        reg('gammamix', ('cpmix', 'R'), lambda u: 'cpmix/(cpmix - R)')
        reg('a', ('gammamix', 'R', 'T'), lambda u: 'sqrt(gammamix*R*T)')
        reg('h', ('T',),
            lambda u: [f'T*{cp0s[n]}' for n in range(ns)])
        reg('e_sp', ('T',),
            lambda u: [f'{sp.cp0 - RU/sp.MW}*T' for sp in self.species])

        # Physical mixture entropy (clamped for non-physical states)
        fpmax = self.fpdtype_max
        cvks = [sp.cp0 - RU/sp.MW for sp in self.species]
        Rks = [RU/sp.MW for sp in self.species]

        def dev_s(u, sfx):
            terms = ' + '.join(
                f'(({u}[{n}] > 0)'
                f' ? Y{sfx}[{n}]*({cvks[n]}*log(T{sfx})'
                f' - {Rks[n]}*log({u}[{n}])) : 0.0)'
                for n in range(ns)
            )

            return (f'fpdtype_t s{sfx} = (T{sfx} > 0)'
                    f' ? ({terms}) : {fpmax};')

        def host_s(u, ns_):
            Y, T = ns_['Y'], ns_['T']
            Tc = np.maximum(T, np.finfo(float).tiny)

            s = sum(np.where(u[n] > 0,
                             Y[n]*(cvks[n]*np.log(Tc)
                                   - Rks[n]*np.log(np.maximum(
                                       u[n], np.finfo(float).tiny))),
                             0.0)
                    for n in range(ns))

            return np.where(T > 0, s, fpmax)

        reg('s', ('T', 'Y'), device=dev_s, host=host_s)

    def pri_to_con(self, pris):
        ndims = self.ndims
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


class MCTPGFluid(MCBaseFluid):
    name = 'mc-tpg'
    species_cls = TPGSpecies

    # Per-species enthalpy at an arbitrary temperature symbol
    def h_of_T(self, n, T):
        return self.species[n].h_expr(T)

    # Newton iteration counts (mirrors mc/develop strict mode)
    dev_niter = 7
    host_niter = 10

    def _eval_cp(self, sp, T):
        cp = eval_over_ranges(T, sp.thermo_ranges, sp.thermo_coeffs,
                              eval_nasa7_cp)
        return cp*RU/sp.MW

    def _eval_h(self, sp, T):
        h = eval_over_ranges(T, sp.thermo_ranges, sp.thermo_coeffs,
                             eval_nasa7_h)
        return h*RU/sp.MW

    def _T0(self):
        T_lo = max(sp.thermo_ranges[0] for sp in self.species)
        T_hi = min(sp.thermo_ranges[-1] for sp in self.species)

        return 0.5*(T_lo + T_hi)

    def _register_thermo(self):
        species = self.species
        reg = self._register_quantity

        def dev_T(u, s):
            lines = [f'fpdtype_t T{s} = {self._T0()};']
            for i in range(self.dev_niter):
                body = ['fpdtype_t h_ = 0.0, cp_ = 0.0;']
                for n, sp in enumerate(species):
                    body.append(f'{{ fpdtype_t cps = {sp.cp_expr(f"T{s}")};')
                    body.append(f'  fpdtype_t hs = {sp.h_expr(f"T{s}")};')
                    body.append(f'  cp_ += cps*Y{s}[{n}];'
                                f' h_ += hs*Y{s}[{n}]; }}')
                body.append(f'T{s} -= (e{s} - (h_ - R{s}*T{s}))'
                            f'/(-cp_ + R{s});')
                lines.append('{\n' + '\n'.join(body) + '\n}')

            return '\n'.join(lines)

        def host_T(u, ns_):
            Y, R, e = ns_['Y'], ns_['R'], ns_['e']

            T = np.zeros_like(np.asarray(e, dtype=float)) + self._T0()
            for _ in range(self.host_niter):
                h = cp = 0.0
                for n, sp in enumerate(species):
                    cp = cp + Y[n]*self._eval_cp(sp, T)
                    h = h + Y[n]*self._eval_h(sp, T)
                T = T - (e - (h - R*T))/(-cp + R)

            return T

        def dev_cpmix(u, s):
            terms = ' + '.join(f'Y{s}[{n}]*({sp.cp_expr(f"T{s}")})'
                               for n, sp in enumerate(species))
            return f'fpdtype_t cpmix{s} = {terms};'

        def host_cpmix(u, ns_):
            Y, T = ns_['Y'], ns_['T']

            return sum(Y[n]*self._eval_cp(sp, T)
                       for n, sp in enumerate(species))

        def dev_h(u, s):
            lines = [f'fpdtype_t h{s}[{self.ns}];']
            lines.extend(f'h{s}[{n}] = {sp.h_expr(f"T{s}")};'
                         for n, sp in enumerate(species))

            return '\n'.join(lines)

        def host_h(u, ns_):
            T = ns_['T']

            return [self._eval_h(sp, T) for sp in species]

        reg('T', ('e', 'R', 'Y'), device=dev_T, host=host_T)
        reg('cpmix', ('T', 'Y'), device=dev_cpmix, host=host_cpmix)
        reg('p', ('rho', 'R', 'T'), lambda u: 'rho*R*T')
        reg('gammamix', ('cpmix', 'R'), lambda u: 'cpmix/(cpmix - R)')
        reg('a', ('gammamix', 'R', 'T'), lambda u: 'sqrt(gammamix*R*T)')
        reg('h', ('T',), device=dev_h, host=host_h)
        reg('e_sp', ('h', 'T'),
            lambda u: [f'h[{n}] - {RU/sp.MW}*T'
                       for n, sp in enumerate(self.species)])

        # Physical mixture entropy (clamped for non-physical states)
        fpmax = self.fpdtype_max
        Rks = [RU/sp.MW for sp in species]

        def dev_s(u, sfx):
            terms = ' + '.join(
                f'(({u}[{n}] > 0)'
                f' ? Y{sfx}[{n}]*(({sp.s_expr(f"T{sfx}")})'
                f' - {Rks[n]}*log({u}[{n}])) : 0.0)'
                for n, sp in enumerate(species)
            )

            return (f'fpdtype_t s{sfx} = (T{sfx} > 0)'
                    f' ? ({terms}) : {fpmax};')

        def host_s(u, ns_):
            Y, T = ns_['Y'], ns_['T']
            Tc = np.maximum(T, np.finfo(float).tiny)

            s = sum(np.where(u[n] > 0,
                             Y[n]*(eval_over_ranges(Tc, sp.thermo_ranges,
                                                    sp.thermo_coeffs,
                                                    eval_nasa7_s)*RU/sp.MW
                                   - Rks[n]*np.log(np.maximum(
                                       u[n], np.finfo(float).tiny))),
                             0.0)
                    for n, sp in enumerate(species))

            return np.where(T > 0, s, fpmax)

        reg('s', ('T', 'Y'), device=dev_s, host=host_s)

    def pri_to_con(self, pris):
        ndims = self.ndims
        p, T = np.asarray(pris[0]), np.asarray(pris[ndims + 1])
        vs = pris[1:ndims + 1]
        Yk = list(pris[ndims + 2:])
        Yns = 1.0 - sum(Yk)
        Yk.append(Yns)

        Rmix = sum(Y / sp.MW for Y, sp in zip(Yk, self.species)) * RU
        h = sum(Y * self._eval_h(sp, T)
                for Y, sp in zip(Yk, self.species))

        rho = p / (Rmix * T)
        rhovs = [rho * v for v in vs]
        rhoE = rho * h - p + 0.5 * rho * sum(v*v for v in vs)
        rhoYk = [rho * Y for Y in Yk]

        return [*rhoYk, *rhovs, rhoE]

    def con_to_pri(self, cons):
        ns = self.ns
        ndims = len(cons) - ns - 1

        rhoY = cons[:ns]
        rho = sum(rhoY)
        vs = [rhov / rho for rhov in cons[ns:ns + ndims]]
        Yk = [rhoYk / rho for rhoYk in rhoY]

        Rmix = sum(Y / sp.MW for Y, sp in zip(Yk, self.species)) * RU

        T = self.eval('T', cons)['T']
        p = rho * Rmix * T

        return [p, *vs, T, *Yk[:-1]]


class MCKineticTheoryTransport(BaseTransport):
    name = 'kinetic-theory'

    def __init__(self, cfg, fluid):
        from pyfr.fluids.collision import astar_at, fit_poly, omega22_at
        from pyfr.fluids.constants import AVOGADRO, EPS0, KB

        if not isinstance(fluid, MCTPGFluid):
            raise ValueError('kinetic-theory transport requires eos = mc-tpg')

        sects = fluid._species_sects
        for name, sect in sects.items():
            if 'transport' not in sect:
                raise ValueError(f'Species {name!r} has no transport data; '
                                 'kinetic-theory transport requires it')

        self.mixing_rule = cfg.get('multi-component', 'mixing-rule', 'Wilke')
        if self.mixing_rule not in ('Wilke', 'Herning-Zipperer'):
            raise ValueError(f'Invalid mixing rule {self.mixing_rule!r}')

        prec = cfg.get('backend', 'precision', 'double')
        self.fpdtype_eps = float(np.finfo(
            np.float32 if prec == 'single' else np.float64).eps)

        tol = cfg.getfloat('multi-component', 'fit-tol', 0.0)

        species = fluid.species
        ns = len(species)
        self.MWs = [sp.MW for sp in species]

        # Per-species Lennard-Jones and internal-dof parameters
        trans = [sects[sp.name]['transport'] for sp in species]
        masses = [sp.MW/AVOGADRO for sp in species]
        diams = [t['diameter'] for t in trans]
        wells = [t['well-depth'] for t in trans]
        dipoles = [t['dipole'] for t in trans]
        polarizs = [t['polarizability'] for t in trans]
        rot_relaxs = [t['rotational-relaxation'] for t in trans]
        rot_dofs = [{'atom': 0.0, 'linear': 1.0}.get(t['geometry'], 1.5)
                    for t in trans]
        polars = [d > 0 for d in dipoles]
        dstar_self = [0.5*d*d/(4*np.pi*EPS0*w*dm**3)
                      for d, w, dm in zip(dipoles, wells, diams)]

        # Global temperature range for the transport fits
        tmin = cfg.get('multi-component', 't-min', '')
        tmax = cfg.get('multi-component', 't-max', '')
        T_lo = (float(tmin) if tmin
                else max(sp.thermo_ranges[0] for sp in species))
        T_hi = (float(tmax) if tmax
                else min(sp.thermo_ranges[-1] for sp in species))

        Ts = np.linspace(T_lo, T_hi, 50)
        logTs = np.log(Ts)
        sqrtTs = np.sqrt(Ts)

        # Pure-species viscosity and conductivity fits
        self.mu_polys, self.kappa_polys = [], []
        for n, sp in enumerate(species):
            o22 = np.array([omega22_at(T*KB/wells[n], dstar_self[n])
                            for T in Ts])
            ast = np.array([astar_at(T*KB/wells[n], dstar_self[n])
                            for T in Ts])
            o11 = o22/ast

            mu = ((5.0/16.0)*np.sqrt(np.pi*masses[n]*KB*Ts)
                  / (np.pi*diams[n]**2*o22))

            r_mass = masses[n]/2.0
            diff_self = ((3.0/16.0)*np.sqrt(2.0*np.pi/r_mass)*(KB*Ts)**1.5
                         / (np.pi*diams[n]**2*o11))

            f_int = sp.MW/(RU*Ts)*diff_self/mu
            Tstar = Ts*KB/wells[n]
            Tstar_298 = KB*298.0/wells[n]
            fz_298 = (1.0 + np.pi**1.5/np.sqrt(Tstar_298)*(0.5
                      + 1.0/Tstar_298) + (0.25*np.pi**2 + 2)/Tstar_298)
            fz_T = (1.0 + np.pi**1.5/np.sqrt(Tstar)*(0.5 + 1.0/Tstar)
                    + (0.25*np.pi**2 + 2)/Tstar)
            A_factor = 2.5 - f_int
            B_factor = (rot_relaxs[n]*fz_298/fz_T
                        + 2.0/np.pi*(5.0/3.0*rot_dofs[n] + f_int))
            c1 = 2.0/np.pi*A_factor/B_factor

            cp_R = eval_over_ranges(Ts, sp.thermo_ranges, sp.thermo_coeffs,
                                    eval_nasa7_cp)
            cv_int = cp_R - 2.5 - rot_dofs[n]
            f_trans = 2.5*(1.0 - c1*rot_dofs[n]/1.5)
            f_rot = f_int*(1.0 + c1)
            kappa = (mu/sp.MW*RU
                     * (f_trans*1.5 + f_rot*rot_dofs[n] + f_int*cv_int))

            mu_scaled = np.sqrt(mu/sqrtTs)
            self.mu_polys.append(fit_poly(logTs, mu_scaled, tol,
                                          w=1.0/np.abs(mu_scaled)))

            kappa_scaled = kappa/sqrtTs
            self.kappa_polys.append(fit_poly(logTs, kappa_scaled, tol,
                                             w=1.0/np.abs(kappa_scaled)))

        # Pairwise binary diffusion fits (of 1/(D/T^1.5) in log T)
        self.dij_polys = [[None]*ns for _ in range(ns)]
        for n in range(ns):
            for m in range(n, ns):
                r_mass = masses[n]*masses[m]/(masses[n] + masses[m])
                r_diam = 0.5*(diams[n] + diams[m])
                r_well = np.sqrt(wells[n]*wells[m])
                r_dipole = np.sqrt(dipoles[n]*dipoles[m])

                r_dstar = 0.5*r_dipole**2/(4*np.pi*EPS0*r_well*r_diam**3)

                if polars[n] != polars[m]:
                    kp, knp = (n, m) if polars[n] else (m, n)
                    alphastar = polarizs[knp]/diams[knp]**3
                    dipolestar = dipoles[kp]/np.sqrt(
                        4*np.pi*EPS0*diams[kp]**3*wells[kp])
                    xi = (1.0 + 0.25*alphastar*dipolestar**2
                          * np.sqrt(wells[kp]/wells[knp]))
                    r_well *= xi**2
                    r_diam *= xi**(-1/6)

                diff = np.empty_like(Ts)
                for i, T in enumerate(Ts):
                    o22 = omega22_at(T*KB/r_well, r_dstar)
                    o11 = o22/astar_at(T*KB/r_well, r_dstar)

                    diff[i] = ((3.0/16.0)*np.sqrt(2.0*np.pi/r_mass)
                               * (KB*T)**1.5/(np.pi*r_diam**2*o11))

                inv_diff_scaled = 1.0/(diff/Ts**1.5)
                poly = fit_poly(logTs, inv_diff_scaled, tol,
                                w=1.0/np.abs(inv_diff_scaled))

                self.dij_polys[n][m] = self.dij_polys[m][n] = poly

    def register(self, fluid):
        from pyfr.fluids.species import BaseSpecies

        ns = len(self.MWs)
        MWs, eps = self.MWs, self.fpdtype_eps
        mu_polys, kappa_polys = self.mu_polys, self.kappa_polys
        dij_polys, mixing_rule = self.dij_polys, self.mixing_rule
        horner = BaseSpecies.horner
        reg = fluid._register_quantity

        # Clamped mole fractions
        def dev_X(u, s):
            lines = [f'fpdtype_t X{s}[{ns}];',
                     '{', 'fpdtype_t total_ = 0.0;']
            for n in range(ns):
                lines.append(f'X{s}[{n}] = Y{s}[{n}]*{1.0/MWs[n]};')
                lines.append(f'total_ += X{s}[{n}];')
            lines.append('fpdtype_t invtotal_ = 1.0/total_;')
            for n in range(ns):
                lines.append(f'X{s}[{n}] = fmax(X{s}[{n}]*invtotal_, {eps});')
            lines.append('}')

            return '\n'.join(lines)

        def host_X(u, ns_):
            Y = ns_['Y']
            Xs = [Y[n]/MWs[n] for n in range(ns)]
            total = sum(Xs)

            return [np.maximum(X/total, eps) for X in Xs]

        reg('X', ('Y',), device=dev_X, host=host_X)
        reg('MWmix', ('X',),
            lambda u: ' + '.join(f'X[{n}]*{MWs[n]}' for n in range(ns)))
        reg('logT', ('T',), lambda u: 'log(T)')
        reg('sqrtT', ('T',), lambda u: 'sqrt(T)')

        # Wilke (temperature-dependent phi) or Herning-Zipperer viscosity
        def dev_mu(u, s):
            lines = [f'fpdtype_t polymu{s}[{ns}], invpolymu{s}[{ns}];']
            for n in range(ns):
                lines.append(f'polymu{s}[{n}] = '
                             f'{horner(f"logT{s}", mu_polys[n])};')
                lines.append(f'invpolymu{s}[{n}] = 1.0/polymu{s}[{n}];')
            lines.append(f'fpdtype_t mu{s} = 0.0;')
            lines.append('{')

            if mixing_rule == 'Wilke':
                lines.append(f'fpdtype_t phi_[{ns}] = {{0}};')
                for n in range(ns):
                    for m in range(ns):
                        if n == m:
                            lines.append(f'phi_[{n}] += X{s}[{m}];')
                        else:
                            fac = 1.0/(np.sqrt(8.0)
                                       * np.sqrt(1.0 + MWs[n]/MWs[m]))
                            lines.append(
                                f'{{ fpdtype_t num_ = 1.0 + polymu{s}[{n}]'
                                f'*invpolymu{s}[{m}]'
                                f'*{(MWs[m]/MWs[n])**0.25};\n'
                                f'phi_[{n}] += num_*num_*{fac}'
                                f'*X{s}[{m}]; }}'
                            )
                for n in range(ns):
                    lines.append(f'mu{s} += polymu{s}[{n}]*polymu{s}[{n}]'
                                 f'*sqrtT{s}*X{s}[{n}]/phi_[{n}];')
            else:
                lines.append('fpdtype_t num_ = 0.0, den_ = 0.0;')
                for n in range(ns):
                    lines.append(
                        f'{{ fpdtype_t xsm_ = X{s}[{n}]'
                        f'*{np.sqrt(MWs[n])};\n'
                        f'num_ += xsm_*polymu{s}[{n}]*polymu{s}[{n}];\n'
                        f'den_ += xsm_; }}'
                    )
                lines.append(f'mu{s} = (num_/den_)*sqrtT{s};')

            lines.append('}')

            return '\n'.join(lines)

        def host_mu(u, ns_):
            X, logT, sqrtT = ns_['X'], ns_['logT'], ns_['sqrtT']
            polymu = [np.polyval(mu_polys[n][::-1], logT) for n in range(ns)]

            if mixing_rule == 'Wilke':
                mu = 0.0
                for n in range(ns):
                    phi = 0.0
                    for m in range(ns):
                        if n == m:
                            phi = phi + X[m]
                        else:
                            fac = 1.0/(np.sqrt(8.0)
                                       * np.sqrt(1.0 + MWs[n]/MWs[m]))
                            num = (1.0 + polymu[n]/polymu[m]
                                   * (MWs[m]/MWs[n])**0.25)
                            phi = phi + num*num*fac*X[m]
                    mu = mu + polymu[n]*polymu[n]*sqrtT*X[n]/phi
            else:
                num = sum(X[n]*np.sqrt(MWs[n])*polymu[n]*polymu[n]
                          for n in range(ns))
                den = sum(X[n]*np.sqrt(MWs[n]) for n in range(ns))
                mu = (num/den)*sqrtT

            return mu

        reg('mu', ('X', 'logT', 'sqrtT'), device=dev_mu, host=host_mu)

        def dev_kappa(u, s):
            lines = [f'fpdtype_t kappa{s};', '{',
                     'fpdtype_t sum1_ = 0.0, sum2_ = 0.0;']
            for n in range(ns):
                lines.append(
                    f'{{ fpdtype_t ks_ = ({horner(f"logT{s}",
                                                  kappa_polys[n])})'
                    f'*sqrtT{s};\n'
                    f'sum1_ += X{s}[{n}]*ks_;\n'
                    f'sum2_ += X{s}[{n}]/ks_; }}'
                )
            lines.append(f'kappa{s} = 0.5*(sum1_ + 1.0/sum2_);')
            lines.append('}')

            return '\n'.join(lines)

        def host_kappa(u, ns_):
            X, logT, sqrtT = ns_['X'], ns_['logT'], ns_['sqrtT']

            sum1 = sum2 = 0.0
            for n in range(ns):
                ks = np.polyval(kappa_polys[n][::-1], logT)*sqrtT
                sum1 = sum1 + X[n]*ks
                sum2 = sum2 + X[n]/ks

            return 0.5*(sum1 + 1.0/sum2)

        reg('kappa', ('X', 'logT', 'sqrtT'), device=dev_kappa,
            host=host_kappa)

        # Mixture-averaged diffusion coefficients
        def dev_D(u, s):
            lines = [f'fpdtype_t D{s}[{ns}];', '{',
                     f'fpdtype_t Tm32_ = 1.0/(sqrtT{s}*sqrtT{s}*sqrtT{s});',
                     f'fpdtype_t sum1_[{ns}] = {{0}};',
                     f'fpdtype_t sum2_[{ns}] = {{0}};']
            for n in range(ns):
                for m in range(n + 1, ns):
                    lines.append(
                        f'{{ fpdtype_t invDij_ = '
                        f'({horner(f"logT{s}", dij_polys[n][m])})*Tm32_;\n'
                        f'fpdtype_t tmp_ = X{s}[{m}]*invDij_;\n'
                        f'sum1_[{n}] += tmp_;\n'
                        f'sum2_[{n}] += tmp_*{MWs[m]};\n'
                        f'tmp_ = X{s}[{n}]*invDij_;\n'
                        f'sum1_[{m}] += tmp_;\n'
                        f'sum2_[{m}] += tmp_*{MWs[n]}; }}'
                    )
            for n in range(ns):
                lines.append(f'sum2_[{n}] *= X{s}[{n}]'
                             f'/(MWmix{s} - {MWs[n]}*X{s}[{n}]);')
                lines.append(f'D{s}[{n}] = '
                             f'1.0/(p{s}*(sum1_[{n}] + sum2_[{n}]));')
            lines.append('}')

            return '\n'.join(lines)

        def host_D(u, ns_):
            X, logT, sqrtT = ns_['X'], ns_['logT'], ns_['sqrtT']
            p, MWmix = ns_['p'], ns_['MWmix']

            Tm32 = 1.0/(sqrtT*sqrtT*sqrtT)
            sum1 = [0.0]*ns
            sum2 = [0.0]*ns
            for n in range(ns):
                for m in range(n + 1, ns):
                    invDij = np.polyval(dij_polys[n][m][::-1], logT)*Tm32
                    sum1[n] = sum1[n] + X[m]*invDij
                    sum2[n] = sum2[n] + X[m]*invDij*MWs[m]
                    sum1[m] = sum1[m] + X[n]*invDij
                    sum2[m] = sum2[m] + X[n]*invDij*MWs[n]

            return [1.0/(p*(sum1[n]
                            + sum2[n]*X[n]/(MWmix - MWs[n]*X[n])))
                    for n in range(ns)]

        reg('D', ('X', 'logT', 'sqrtT', 'p', 'MWmix'), device=dev_D,
            host=host_D)


class MCConstantPropsTransport(BaseTransport):
    name = 'constant-props'

    def __init__(self, cfg, fluid):
        sects = getattr(fluid, '_species_sects', None)
        if sects is None:
            raise ValueError('constant-props transport requires a '
                             'multi-component fluid')

        try:
            self.mu0s = [float(s['mu0']) for s in sects.values()]
            self.kappa0s = [float(s['kappa0']) for s in sects.values()]
        except KeyError as e:
            raise ValueError(f'constant-props transport requires {e} for '
                             'every species (yaml or ini)') from None

        self.Les = [float(s.get('Le', 1.0)) for s in sects.values()]
        self.MWs = [sp.MW for sp in fluid.species]

    def register(self, fluid):
        ns = len(self.MWs)
        eps = sys.float_info.epsilon
        mu0s, kappa0s, MWs, Les = self.mu0s, self.kappa0s, self.MWs, self.Les
        reg = fluid._register_quantity

        # Mole fractions
        denom = ' + '.join(f'Y[{m}]*{1.0 / MWs[m]}' for m in range(ns))
        reg('X', ('Y',),
            lambda u: [f'Y[{n}]*{1.0 / MWs[n]}/({denom})' for n in range(ns)])

        # Wilke mixture viscosity; the phi coefficients are constants
        phi = [[(1 + math.sqrt(mu0s[n] / max(mu0s[m], eps)
                               * math.sqrt(max(mu0s[m], eps) / mu0s[n])))**2
                / (math.sqrt(8)*math.sqrt(1 + MWs[n]/MWs[m]))
                for m in range(ns)] for n in range(ns)]

        mu_terms = ' + '.join(
            f'{mu0s[n]}*X[{n}]/('
            + ' + '.join(f'{phi[n][m]}*X[{m}]' for m in range(ns)) + ')'
            for n in range(ns)
        )
        reg('mu', ('X',), lambda u: mu_terms)

        # Mixture thermal conductivity
        sum1 = ' + '.join(f'X[{n}]*{kappa0s[n]}' for n in range(ns))
        sum2 = ' + '.join(f'X[{n}]*{1.0 / (eps + kappa0s[n])}'
                          for n in range(ns))
        reg('kappa', ('X',), lambda u: f'0.5*(({sum1}) + 1.0/({sum2}))')

        # Lewis number approximation for the species diffusivities
        reg('D', ('kappa', 'rho', 'cpmix'),
            lambda u: [f'kappa/(rho*cpmix*{Les[n]})' for n in range(ns)])
