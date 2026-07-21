import numpy as np

from pyfr.multicomp import RU
from pyfr.multicomp.chem.reaction import Reaction
from pyfr.util import subclass_where


class Chemistry:
    def __init__(self, cfg, *args, **kwargs):
        super().__init__(cfg, *args, **kwargs)

        self.reactions = []
        for i, rxn_sect in enumerate(self._reaction_sects):
            rtype = rxn_sect['rtype']
            rxn_cls = subclass_where(Reaction, name=rtype)
            self.reactions.append(rxn_cls(i, rxn_sect))

        # Populate each reaction with species data for self-contained expressions
        nu_f = self.nu_f
        nu_b = self.nu_b
        aij = self.aij
        MWs = self.MWs
        for rxn in self.reactions:
            rxn.populate(self.ns, MWs, nu_f[rxn.index], nu_b[rxn.index],
                         aij[rxn.index], self.sp_idx)

    @property
    def nr(self):
        return len(self.reactions)

    def reactions_by_type(self, name):
        return [rxn for rxn in self.reactions if rxn.name == name]

    @property
    def has_reversible(self):
        return any(r.reversible and not r.zero_rate for r in self.reactions)

    @property
    def has_third_body(self):
        return any(hasattr(r, 'efficiencies') and not r.zero_rate
                   for r in self.reactions)

    @property
    def gbs_species(self):
        # Species whose Gibbs energy enters some reversible reaction
        s = set()
        for r in self.reactions:
            if r.reversible and not r.zero_rate:
                s.update(n for n, v in enumerate(r.nu_sum) if v != 0)
        return sorted(s)

    @property
    def conc_species(self):
        # Species whose log-concentration is referenced by any rate
        # expression (forward orders/stoich, or reverse net stoich)
        s = set()
        for r in self.reactions:
            if r.zero_rate:
                continue
            s.update(r.fwd_species)
            if r.reversible:
                s.update(n for n, v in enumerate(r.nu_sum) if v != 0)
        return sorted(s)

    @property
    def ctot_coeffs(self):
        # Inverse molecular weights for the shared total concentration
        return [(n, 1.0/mw) for n, mw in enumerate(self.MWs)]

    @property
    def max_omega_gain(self):
        # Worst-case amplification from a single rate of progress to a
        # mass production rate: max_n MW_n * sum_i |nu_in|
        gain = np.zeros(self.ns)
        for r in self.reactions:
            if not r.zero_rate:
                for n, v in r.nu_terms:
                    gain[n] += abs(v)
        return float(max(1.0, (gain*self.MWs).max()))

    @property
    def log_c_floor(self):
        """Log-concentration floor for absent species, derived from
        the mechanism.

        A species with zero concentration must suppress each direction
        of every reaction it participates in: for a floored species
        with exponent v in a direction whose remaining factors can
        reach a log-rate of B, we need v*|floor| >= B + margin.  The
        temperature-dependent parts (ln kf, the equilibrium affinity)
        are evaluated over the mechanism's own thermo fit range; the
        remaining leaf bounds are ln c <= 12 per concentration factor
        (c <= 1.6e5 kmol/m^3, beyond any gas state) and a margin
        pushing suppressed rates below e^-45.  The floor is kept as
        shallow as the mechanism allows because low-precision builds
        round the rate exponent at one ulp of its largest intermediate.
        """
        log_c_max = 12.0
        margin = 45.0

        tlo = min((sp.thermo_ranges[0] for sp in self.species
                   if hasattr(sp, 'thermo_ranges')), default=200.0)
        thi = max((sp.thermo_ranges[-1] for sp in self.species
                   if hasattr(sp, 'thermo_ranges')), default=5000.0)
        T = np.geomspace(max(tlo, 100.0), thi, 256)
        logT = np.log(T)

        gbs = {}

        def gbs_of(n):
            if n not in gbs:
                gbs[n] = self.species[n].gbs_eval(T)
            return gbs[n]

        # Never shallower than the legacy FLT_MIN concentration floor
        req = 90.0
        for r in self.reactions:
            if r.zero_rate:
                continue
            lkf = r.rate.log_A + r.rate.b*logT - r.rate.Ea/T
            # Three-body rates carry one cTBC factor; falloff F <= ~2
            if r.name == 'three-body':
                aux = log_c_max + 2.0
            elif r.name.startswith('falloff'):
                aux = 1.0
            else:
                aux = 0.0

            fmax_ = lkf.max() + aux
            exps = r.fwd_exps
            for k, vk in exps:
                if vk <= 0:
                    continue  # negative orders cannot be suppressed
                b = fmax_ + sum(v*log_c_max for n, v in exps
                                if n != k and v > 0)
                req = max(req, (b + margin)/vk)

            if r.reversible:
                aff = sum(v*gbs_of(n) for n, v in r.nu_terms)
                lkr = lkf + aff + r.nu_total*np.log(RU*T/101325.0)
                rexps = [(n, float(v)) for n, v in enumerate(r.nu_b_row)
                         if v > 0]
                rmax = lkr.max() + aux
                for k, vk in rexps:
                    b = rmax + sum(v*log_c_max for n, v in rexps if n != k)
                    req = max(req, (b + margin)/vk)

        return -float(req)

    @property
    def nu_f(self):
        mat = np.zeros((self.nr, self.ns))
        for rxn in self.reactions:
            for name, coeff in rxn.reactants.items():
                mat[rxn.index, self.sp_idx(name)] = coeff
        return mat

    @property
    def nu_b(self):
        mat = np.zeros((self.nr, self.ns))
        for rxn in self.reactions:
            for name, coeff in rxn.products.items():
                mat[rxn.index, self.sp_idx(name)] = coeff
        return mat

    @property
    def aij(self):
        mat = np.ones((self.nr, self.ns))
        for rxn in self.reactions:
            if hasattr(rxn, 'efficiencies'):
                mat[rxn.index, :] = getattr(rxn, 'default_efficiency', 1.0)
                for name, eff in rxn.efficiencies.items():
                    mat[rxn.index, self.sp_idx(name)] = eff
        return mat
