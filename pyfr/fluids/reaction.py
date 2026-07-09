import math

from pyfr.fluids.species import clean_csigns


class ArrheniusRate:
    def __init__(self, rate_sect):
        self.A = float(rate_sect['A'])
        self.b = float(rate_sect['b'])
        self.Ea = float(rate_sect['Ea'])


class Reaction:
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

    def populate(self, ns, MWs, nu_f_row, nu_b_row, aij_row, sp_idx_fn):
        self.ns = ns
        self.MWs = MWs
        self.nu_f_row = nu_f_row
        self.nu_b_row = nu_b_row
        self.nu_sum = nu_b_row - nu_f_row
        self.aij_row = aij_row
        if self.orders:
            self.orders = {sp_idx_fn(k): v for k, v in self.orders.items()}

    def log_kf_expr(self, vlogT='logT', vTinv='Tinv'):
        r = self.rate
        return f'{math.log(r.A)} + ({r.b}*{vlogT}) - ({r.Ea}*{vTinv})'

    def fwd_rate_expr(self, vlog_cs='log_cs'):
        if self.orders:
            terms = [f'({v})*{vlog_cs}[{n}]'
                     for n, v in self.orders.items() if float(v) != 0]
        else:
            terms = [f'({v})*{vlog_cs}[{n}]'
                     for n, v in enumerate(self.nu_f_row) if v != 0]
        return ' + '.join(terms) if terms else '0'

    def reverse_rate_block(self, vlog_cs='log_cs', vgbs='gbs',
                           vprefRuT='log_prefRuT',
                           vprefRuTinv='log_prefRuTinv'):
        lines = []
        # Equilibrium constant from Gibbs
        kp_terms = [f'({v})*{vgbs}[{n}]'
                    for n, v in enumerate(self.nu_sum) if v != 0]
        lines.append(f'fpdtype_t log_Kp = {" + ".join(kp_terms)};')

        # Kc correction for concentration units
        nu_total = float(sum(self.nu_sum))
        if nu_total > 0:
            log_term = f'{nu_total}*{vprefRuTinv}'
        elif nu_total < 0:
            log_term = f'{-nu_total}*{vprefRuT}'
        else:
            log_term = '0.0'
        lines.append(f'fpdtype_t log_k_r = log_k_f + log_Kp + {log_term};')

        # Reverse rate: use expm1 to avoid cancellation near equilibrium
        rb_terms = [f'({v})*{vlog_cs}[{n}]'
                    for n, v in enumerate(self.nu_b_row) if v != 0]
        lines.append(f'fpdtype_t log_rp_reverse = log_k_r + {" + ".join(rb_terms)};')
        lines.append('fpdtype_t log_rp_max = fmax(log_rp, log_rp_reverse);')
        lines.append('fpdtype_t rp_diff = log_rp - log_rp_reverse;')
        lines.append('rp = copysign(exp(log_rp_max), rp_diff) * (-expm1(-fabs(rp_diff)));')
        return '\n    '.join(lines)

    def accumulate_omega_expr(self, vomega='omega'):
        lines = []
        for n, nu in enumerate(self.nu_sum):
            if abs(nu) > 0:
                lines.append(f'{vomega}[{n}] += {nu}*rp;')
        return '\n    '.join(lines)


class ElementaryReaction(Reaction):
    name = 'elementary'

    @clean_csigns
    def rate_block(self, vlogT='logT', vTinv='Tinv', vlog_cs='log_cs',
                   vgbs='gbs', vrho='rho', vY='Y', vomega='omega',
                   vprefRuT='log_prefRuT', vprefRuTinv='log_prefRuTinv'):
        lines = [f'// R{self.index}: {self.equation}']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr(vlogT, vTinv)};')
        lines.append(f'  fpdtype_t log_rp = log_k_f + {self.fwd_rate_expr(vlog_cs)};')
        lines.append('  fpdtype_t rp = exp(log_rp);')
        if self.reversible:
            lines.append(f'  {self.reverse_rate_block(vlog_cs, vgbs, vprefRuT, vprefRuTinv)}')
        lines.append(f'  {self.accumulate_omega_expr(vomega)}')
        lines.append('}')
        return '\n'.join(lines)


class ThreeBodyReaction(Reaction):
    name = 'three-body'

    def __init__(self, index, rxn_sect):
        super().__init__(index, rxn_sect)
        self.efficiencies = dict(rxn_sect.get('efficiencies', {}))

    def ctbc_expr(self, vrho='rho', vY='Y'):
        terms = [f'({eff/self.MWs[n]})*{vY}[{n}]'
                 for n, eff in enumerate(self.aij_row) if eff != 0]
        return f'{vrho} * ({"+".join(terms)})'

    @clean_csigns
    def rate_block(self, vlogT='logT', vTinv='Tinv', vlog_cs='log_cs',
                   vgbs='gbs', vrho='rho', vY='Y', vomega='omega',
                   vprefRuT='log_prefRuT', vprefRuTinv='log_prefRuTinv'):
        lines = [f'// R{self.index}: {self.equation}']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr(vlogT, vTinv)};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr(vrho, vY)};')
        lines.append('  log_k_f += log(cTBC);')
        lines.append(f'  fpdtype_t log_rp = log_k_f + {self.fwd_rate_expr(vlog_cs)};')
        lines.append('  fpdtype_t rp = exp(log_rp);')
        if self.reversible:
            lines.append(f'  {self.reverse_rate_block(vlog_cs, vgbs, vprefRuT, vprefRuTinv)}')
        lines.append(f'  {self.accumulate_omega_expr(vomega)}')
        lines.append('}')
        return '\n'.join(lines)


class FalloffReaction(ThreeBodyReaction):
    name = 'falloff-lindemann'

    def __init__(self, index, rxn_sect):
        super().__init__(index, rxn_sect)
        self.low_rate = ArrheniusRate(rxn_sect['low_rate'])

    def log_pr_expr(self, vlogT='logT', vTinv='Tinv'):
        r, lo = self.rate, self.low_rate
        return (f'{math.log(lo.A/r.A)} + '
                f'({lo.b - r.b}*{vlogT}) - '
                f'({lo.Ea - r.Ea}*{vTinv})')

    @clean_csigns
    def rate_block(self, vlogT='logT', vTinv='Tinv', vlog_cs='log_cs',
                   vgbs='gbs', vrho='rho', vY='Y', vomega='omega',
                   vprefRuT='log_prefRuT', vprefRuTinv='log_prefRuTinv'):
        lines = [f'// R{self.index}: {self.equation} (Lindemann)']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr(vlogT, vTinv)};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr(vrho, vY)};')
        lines.append(f'  fpdtype_t log_Pr = log(cTBC) + {self.log_pr_expr(vlogT, vTinv)};')
        lines.append('  fpdtype_t log_pmod = -log1p(exp(-log_Pr));')
        lines.append('  log_k_f += log_pmod;')
        lines.append(f'  fpdtype_t log_rp = log_k_f + {self.fwd_rate_expr(vlog_cs)};')
        lines.append('  fpdtype_t rp = exp(log_rp);')
        if self.reversible:
            lines.append(f'  {self.reverse_rate_block(vlog_cs, vgbs, vprefRuT, vprefRuTinv)}')
        lines.append(f'  {self.accumulate_omega_expr(vomega)}')
        lines.append('}')
        return '\n'.join(lines)


class TroeReaction(FalloffReaction):
    name = 'falloff-troe'

    def __init__(self, index, rxn_sect):
        super().__init__(index, rxn_sect)
        troe = rxn_sect['troe']
        self.troe_A = float(troe['A'])
        self.troe_T3 = float(troe['T3'])
        self.troe_T1 = float(troe['T1'])
        self.troe_T2 = float(troe['T2']) if 'T2' in troe else None

    def _log10_fcent_block(self, vT='T', vTinv='Tinv'):
        alpha = self.troe_A
        Tsss = self.troe_T3
        Ts = self.troe_T1
        Tss = self.troe_T2 if self.troe_T2 is not None else 0.0
        ln10 = 1.0 / math.log(10.0)
        is_3p = Tss == 0.0

        def lse2(a, b, diff):
            op = '-' if diff else '+'
            return (f'{{ fpdtype_t _ref = fmax({a}, {b}); '
                    f'log10Fcent = (_ref + log(exp({a} - _ref) {op} exp({b} - _ref))) * {ln10}; }}')

        def lse3(a, b, c):
            return (f'{{ fpdtype_t _ref = fmax(fmax({a}, {b}), {c}); '
                    f'log10Fcent = (_ref + log(exp({a} - _ref) + exp({b} - _ref) + exp({c} - _ref))) * {ln10}; }}')

        t1 = f'-{vT}*{1.0/Tsss}'
        t2 = f'-{vT}*{1.0/Ts}'
        t3 = f'-{Tss}*{vTinv}' if not is_3p else None

        if alpha == 0.0:
            if is_3p:
                return f'fpdtype_t log10Fcent = -{vT}*{ln10/Tsss};'
            else:
                return f'fpdtype_t log10Fcent;\n  {lse2(t1, t3, False)}'
        elif alpha == 1.0:
            if is_3p:
                return f'fpdtype_t log10Fcent = -{vT}*{ln10/Ts};'
            else:
                return f'fpdtype_t log10Fcent;\n  {lse2(t2, t3, False)}'
        elif alpha > 1.0:
            la = f'{math.log(alpha)} + {t2}'
            lam1 = f'{math.log(alpha - 1.0)} + {t1}'
            if is_3p:
                return f'fpdtype_t log10Fcent;\n  {lse2(la, lam1, True)}'
            else:
                lines = ['fpdtype_t log10Fcent;']
                lines.append(f'  {{ fpdtype_t _ref = fmax({la}, {t3}); '
                             f'fpdtype_t _pos = _ref + log(exp({la} - _ref) + exp({t3} - _ref)); '
                             f'_ref = fmax(_pos, {lam1}); '
                             f'log10Fcent = (_ref + log(exp(_pos - _ref) - exp({lam1} - _ref))) * {ln10}; }}')
                return '\n  '.join(lines)
        elif alpha > 0.0:
            l1ma = f'{math.log(1.0 - alpha)} + {t1}'
            la = f'{math.log(alpha)} + {t2}'
            if is_3p:
                return f'fpdtype_t log10Fcent;\n  {lse2(l1ma, la, False)}'
            else:
                return f'fpdtype_t log10Fcent;\n  {lse3(l1ma, la, t3)}'
        else:  # alpha < 0
            l1ma = f'{math.log(1.0 - alpha)} + {t1}'
            lna = f'{math.log(-alpha)} + {t2}'
            if is_3p:
                return f'fpdtype_t log10Fcent;\n  {lse2(l1ma, lna, True)}'
            else:
                lines = ['fpdtype_t log10Fcent;']
                lines.append(f'  {{ fpdtype_t _ref = fmax({l1ma}, {t3}); '
                             f'fpdtype_t _pos = _ref + log(exp({l1ma} - _ref) + exp({t3} - _ref)); '
                             f'_ref = fmax(_pos, {lna}); '
                             f'log10Fcent = (_ref + log(exp(_pos - _ref) - exp({lna} - _ref))) * {ln10}; }}')
                return '\n  '.join(lines)

    @clean_csigns
    def rate_block(self, vlogT='logT', vTinv='Tinv', vlog_cs='log_cs',
                   vgbs='gbs', vrho='rho', vY='Y', vomega='omega',
                   vprefRuT='log_prefRuT', vprefRuTinv='log_prefRuTinv',
                   vT='T'):
        ln10 = math.log(10.0)
        lines = [f'// R{self.index}: {self.equation} (Troe)']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr(vlogT, vTinv)};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr(vrho, vY)};')
        lines.append(f'  fpdtype_t log_Pr = log(cTBC) + {self.log_pr_expr(vlogT, vTinv)};')
        lines.append(f'  {self._log10_fcent_block(vT, vTinv)}')
        lines.append(f'  fpdtype_t C = -0.4 - 0.67*log10Fcent;')
        lines.append(f'  fpdtype_t N = 0.75 - 1.27*log10Fcent;')
        lines.append(f'  fpdtype_t log10_Pr = log_Pr * {1.0/ln10};')
        lines.append(f'  fpdtype_t A_troe = log10_Pr + C;')
        lines.append(f'  fpdtype_t f1 = A_troe/(N - 0.14*A_troe);')
        lines.append(f'  fpdtype_t log_F_pdr = log10Fcent/(1.0+f1*f1) * {ln10};')
        lines.append('  fpdtype_t log_pmod = -log1p(exp(-log_Pr)) + log_F_pdr;')
        lines.append('  log_k_f += log_pmod;')
        lines.append(f'  fpdtype_t log_rp = log_k_f + {self.fwd_rate_expr(vlog_cs)};')
        lines.append('  fpdtype_t rp = exp(log_rp);')
        if self.reversible:
            lines.append(f'  {self.reverse_rate_block(vlog_cs, vgbs, vprefRuT, vprefRuTinv)}')
        lines.append(f'  {self.accumulate_omega_expr(vomega)}')
        lines.append('}')
        return '\n'.join(lines)
