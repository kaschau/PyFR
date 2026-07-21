import math


class ArrheniusRate:
    def __init__(self, rate_sect):
        self.A = float(rate_sect['A'])
        self.b = float(rate_sect['b'])
        self.Ea = float(rate_sect['Ea'])

    @property
    def log_A(self):
        return math.log(self.A)


class Reaction:
    """Per-reaction mechanism analysis.

    Instances expose plain data (folded constants, exponent and
    stoichiometry tables) consumed by the net-rate-of-production
    template, which owns all of the generated source.
    """

    name = None

    def __init__(self, index, rxn_sect):
        self.index = int(index)
        self.equation = rxn_sect['equation']
        self.reversible = rxn_sect.get('reversible', True)
        self.duplicate = rxn_sect.get('duplicate', False)
        self.reactants = dict(rxn_sect['reactants'])
        self.products = dict(rxn_sect['products'])
        self.rate = ArrheniusRate(rxn_sect['rate'])
        self.orders = dict(rxn_sect.get('orders', {}))

        self._check_rate(self.rate)

    def _check_rate(self, rate):
        # Log-space evaluation requires A > 0; A = 0 zeroes the
        # reaction, which is handled by skipping it entirely
        if rate.A < 0:
            raise ValueError(
                f"R{self.index} '{self.equation}': negative pre-exponential "
                f"factors (negative-A) are not supported by log-space "
                f"chemistry")

    @property
    def zero_rate(self):
        return self.rate.A == 0

    def populate(self, ns, MWs, nu_f_row, nu_b_row, aij_row, sp_idx_fn):
        self.ns = ns
        self.MWs = MWs
        self.nu_f_row = nu_f_row
        self.nu_b_row = nu_b_row
        self.nu_sum = nu_b_row - nu_f_row
        self.aij_row = aij_row
        if self.orders:
            self.orders = {sp_idx_fn(k): v for k, v in self.orders.items()}

    @property
    def fwd_exps(self):
        """Forward concentration exponents as (species, exponent).

        Explicit orders override the stoichiometric order per species;
        unlisted reactants keep their stoichiometric order and
        nonreactant orders enter as additional factors (Cantera rules).
        """
        if self.orders:
            orders = dict(self.orders)
            exps = [(n, orders.pop(n, v))
                    for n, v in enumerate(self.nu_f_row) if v != 0]
            exps += sorted(orders.items())
        else:
            exps = [(n, v) for n, v in enumerate(self.nu_f_row) if v != 0]
        return [(n, float(v)) for n, v in exps if float(v) != 0]

    @property
    def fwd_species(self):
        return {n for n, v in self.fwd_exps}

    @property
    def nu_terms(self):
        """Net stoichiometry as (species, nu_b - nu_f), zeros dropped."""
        return [(n, float(v)) for n, v in enumerate(self.nu_sum) if v != 0]

    @property
    def rev_exps(self):
        """Reverse concentration exponents as (species, nu_b)."""
        return [(n, float(v)) for n, v in enumerate(self.nu_b_row) if v != 0]

    @property
    def nu_total(self):
        return float(sum(self.nu_sum))


class ElementaryReaction(Reaction):
    name = 'elementary'


class ThreeBodyReaction(Reaction):
    name = 'three-body'

    def __init__(self, index, rxn_sect):
        super().__init__(index, rxn_sect)
        self.efficiencies = dict(rxn_sect.get('efficiencies', {}))
        self.default_efficiency = float(
            rxn_sect.get('default-efficiency', 1.0))

    @property
    def tbc_devs(self):
        """Third-body deviation coefficients as (species, (eff - d)/MW).

        Cantera form: concm = default*ctot + sum_k (eff_k - default)*c_k
        with only non-default efficiencies contributing to the sum.
        """
        d = self.default_efficiency
        return [(n, (e - d)/self.MWs[n])
                for n, e in enumerate(self.aij_row) if e != d]


class FalloffReaction(ThreeBodyReaction):
    name = 'falloff-lindemann'

    def __init__(self, index, rxn_sect):
        super().__init__(index, rxn_sect)
        self.low_rate = ArrheniusRate(rxn_sect['low_rate'])
        self._check_rate(self.low_rate)

    @property
    def zero_rate(self):
        return self.rate.A == 0 or self.low_rate.A == 0

    @property
    def pr_rate(self):
        """k0/kinf as a pseudo-Arrhenius rate for log_Pr evaluation."""
        return ArrheniusRate({'A': self.low_rate.A/self.rate.A,
                              'b': self.low_rate.b - self.rate.b,
                              'Ea': self.low_rate.Ea - self.rate.Ea})


class TroeReaction(FalloffReaction):
    name = 'falloff-troe'

    def __init__(self, index, rxn_sect):
        super().__init__(index, rxn_sect)
        troe = rxn_sect['troe']
        self.troe_A = float(troe['A'])
        self.troe_T3 = float(troe['T3'])
        self.troe_T1 = float(troe['T1'])
        self.troe_T2 = float(troe['T2']) if 'T2' in troe else None

    @property
    def fcent_terms(self):
        """Fcent = (1 - a)exp(-T/T3) + a exp(-T/T1) [+ exp(-T2/T)] as
        (positive, const, cT, cTinv) tuples with each term equal to
        +/- exp(const + cT*T + cTinv/T).

        Cantera zeroes the T3/T1 terms below 1e-300 in magnitude —
        indistinguishable from == 0 for any real mechanism — and omits
        the T2 term when it is absent or exactly zero.  Mechanisms
        also disable terms with sentinel values (T3 = 1e-15,
        T1 = 1e+50, ...) whose reciprocals overflow or underflow fp32
        literals; any term whose exponential is identically zero (or
        whose exponent is identically negligible) over gas-phase
        temperatures [100, 20000] K is therefore resolved here.
        """
        tmin, tmax = 100.0, 20000.0
        alpha = self.troe_A
        cand = []
        if alpha != 1.0 and self.troe_T3 != 0:
            c = 1.0 - alpha
            cand.append((c > 0, math.log(abs(c)), -1.0/self.troe_T3, 0.0))
        if alpha != 0.0 and self.troe_T1 != 0:
            cand.append((alpha > 0, math.log(abs(alpha)),
                         -1.0/self.troe_T1, 0.0))
        if self.troe_T2:
            cand.append((True, 0.0, 0.0, -self.troe_T2))

        terms = []
        for pos, c0, cT, cTinv in cand:
            # Decaying exponential that is zero at every temperature
            if cT < 0 and -cT*tmin > 745 or cTinv < 0 and -cTinv/tmax > 745:
                continue
            # Exponent contributions below fp64 resolution
            if cT != 0 and abs(cT)*tmax < 1e-14:
                cT = 0.0
            if cTinv != 0 and abs(cTinv)/tmin < 1e-14:
                cTinv = 0.0
            terms.append((pos, c0, cT, cTinv))
        return terms
