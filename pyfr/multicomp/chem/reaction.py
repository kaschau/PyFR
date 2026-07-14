import math

from pyfr.multicomp import RU, clean_csigns


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
            # Cantera rejects explicit orders on reversible reactions as
            # they break detailed balance; so do we
            if self.reversible:
                raise ValueError('Explicit orders on reversible reaction '
                                 f'{self.equation!r} are not supported')
            self.orders = {sp_idx_fn(k): v for k, v in self.orders.items()}

    def log_kf_expr(self, vlogT='logT', vTinv='Tinv'):
        r = self.rate
        return f'{math.log(r.A)} + ({r.b}*{vlogT}) - ({r.Ea}*{vTinv})'

    def fwd_orders(self):
        # Cantera semantics: explicit orders override the listed species
        # only; unlisted reactants keep their stoichiometric coefficient
        orders = {n: float(v) for n, v in enumerate(self.nu_f_row) if v != 0}
        if self.orders:
            orders |= {n: float(v) for n, v in self.orders.items()}
        return {n: v for n, v in orders.items() if v != 0}

    def fwd_rate_expr(self, vlog_cs='log_cs'):
        terms = [f'({v})*{vlog_cs}[{n}]'
                 for n, v in sorted(self.fwd_orders().items())]
        return ' + '.join(terms) if terms else '0'

    def _reverse_kr_lines(self, vlog_cs='log_cs', vgbs='gbs',
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

        rb_terms = [f'({v})*{vlog_cs}[{n}]'
                    for n, v in enumerate(self.nu_b_row) if v != 0]
        lines.append(f'fpdtype_t log_rp_reverse = log_k_r + {" + ".join(rb_terms)};')
        return lines

    def reverse_rate_block(self, vlog_cs='log_cs', vgbs='gbs',
                           vprefRuT='log_prefRuT',
                           vprefRuTinv='log_prefRuTinv'):
        lines = self._reverse_kr_lines(vlog_cs, vgbs, vprefRuT, vprefRuTinv)

        # Reverse rate: use expm1 to avoid cancellation near equilibrium
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

    def dlogkf_dT_expr(self, vTinv='Tinv'):
        r = self.rate
        return f'({r.b} + {r.Ea}*{vTinv})*{vTinv}'

    def dlogkc_dT_expr(self, vdgbs='dgbs', vTinv='Tinv'):
        # d/dT of (log_Kp + concentration correction) from reverse_rate_block
        terms = [f'({v})*{vdgbs}[{n}]'
                 for n, v in enumerate(self.nu_sum) if v != 0]
        if (nu_total := float(sum(self.nu_sum))) != 0:
            terms.append(f'{nu_total}*{vTinv}')
        return ' + '.join(terms) if terms else '0.0'

    def _jac_rates_lines(self, need_rp=False, vlog_cs='log_cs', vgbs='gbs',
                         vprefRuT='log_prefRuT', vprefRuTinv='log_prefRuTinv'):
        # Forward/reverse rates of progress, kept separate for derivatives;
        # the cancellation-safe net rate rp is only assembled when needed
        # (third-body/falloff common factors)
        lines = [f'fpdtype_t log_rp = log_k_f + {self.fwd_rate_expr(vlog_cs)};']
        lines.append('fpdtype_t q_f = exp(log_rp);')
        if need_rp:
            lines.append('fpdtype_t rp = q_f;')
        if self.reversible:
            if need_rp:
                lines.append(self.reverse_rate_block(vlog_cs, vgbs, vprefRuT,
                                                     vprefRuTinv))
            else:
                lines += self._reverse_kr_lines(vlog_cs, vgbs, vprefRuT,
                                                vprefRuTinv)
            lines.append('fpdtype_t q_r = exp(log_rp_reverse);')
        return lines

    def _jac_dT_lines(self, nvars, extra_dT=None, vdgbs='dgbs', vTinv='Tinv'):
        # d(rp)/dT at fixed concentrations; staged in the energy column
        dlk = self.dlogkf_dT_expr(vTinv)
        if extra_dT:
            dlk = f'{dlk} + {extra_dT}'
        lines = [f'fpdtype_t dlk_dT = {dlk};']
        if self.reversible:
            lines.append(f'fpdtype_t drp_dT = q_f*dlk_dT '
                         f'- q_r*(dlk_dT + {self.dlogkc_dT_expr(vdgbs, vTinv)});')
        else:
            lines.append('fpdtype_t drp_dT = q_f*dlk_dT;')
        for n, nu in enumerate(self.nu_sum):
            if nu != 0:
                lines.append(f'jac[{n*nvars + nvars - 1}] += {nu}*drp_dT;')
        return lines

    def _jac_dC_lines(self, nvars, vsmall, tb_factor=None, vlog_cs='log_cs'):
        # d(rp)/dC_j: kinetic-order terms for participants, plus an optional
        # rank-one third-body/falloff term over all species.  Where the
        # concentration floor is active the coded rate is constant in u_j,
        # so the derivative of the clamped function is zero (this also
        # prevents C^(o-1) blow-up for fractional orders o < 1)
        log_small = math.log(float(vsmall))
        lines = []
        fwd = self.fwd_orders()
        rev = ({n: float(v) for n, v in enumerate(self.nu_b_row) if v != 0}
               if self.reversible else {})
        for j in sorted(set(fwd) | set(rev)):
            terms = []
            if j in fwd:
                terms.append(f'{fwd[j]}*exp(log_rp - {vlog_cs}[{j}])')
            if j in rev:
                terms.append(f'- {rev[j]}*exp(log_rp_reverse - {vlog_cs}[{j}])')
            lines.append(f'{{ fpdtype_t drp_dC = ({vlog_cs}[{j}] > {log_small})'
                         f' ? {" ".join(terms)} : 0.0;')
            for n, nu in enumerate(self.nu_sum):
                if nu != 0:
                    lines.append(f'  jac[{n*nvars + j}] += {nu}*drp_dC;')
            lines.append('}')
        if tb_factor is not None:
            aeff = ', '.join(str(float(a)) for a in self.aij_row)
            lines.append(f'{{ fpdtype_t _aj[{self.ns}] = {{ {aeff} }};')
            lines.append(f'  fpdtype_t _rc = rp*({tb_factor});')
            lines.append(f'  for (int _j = 0; _j < {self.ns}; _j += 1)')
            lines.append('  {')
            for n, nu in enumerate(self.nu_sum):
                if nu != 0:
                    lines.append(f'    jac[{n}*{nvars} + _j] += {nu}*_rc*_aj[_j];')
            lines.append('  }')
            lines.append('}')
        return lines


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


    @clean_csigns
    def jac_block(self, nvars, vsmall='1e-300'):
        lines = [f'// R{self.index}: {self.equation} (Jacobian)']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr()};')
        lines += [f'  {l}' for l in self._jac_rates_lines()]
        lines += [f'  {l}' for l in self._jac_dT_lines(nvars)]
        lines += [f'  {l}' for l in self._jac_dC_lines(nvars, vsmall)]
        lines.append('}')
        return '\n'.join(lines)


class ThreeBodyReaction(Reaction):
    name = 'three-body'

    def __init__(self, index, rxn_sect):
        super().__init__(index, rxn_sect)
        self.efficiencies = dict(rxn_sect.get('efficiencies', {}))
        self.default_efficiency = float(
            rxn_sect.get('default-efficiency', 1.0))

    def ctbc_expr(self, vrho='rho', vY='Y'):
        terms = [f'({eff/self.MWs[n]})*{vY}[{n}]'
                 for n, eff in enumerate(self.aij_row) if eff != 0]
        return f'{vrho} * ({"+".join(terms)})'

    @clean_csigns
    def rate_block(self, vlogT='logT', vTinv='Tinv', vlog_cs='log_cs',
                   vgbs='gbs', vrho='rho', vY='Y', vomega='omega',
                   vprefRuT='log_prefRuT', vprefRuTinv='log_prefRuTinv',
                   vsmall='1e-300'):
        lines = [f'// R{self.index}: {self.equation}']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr(vlogT, vTinv)};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr(vrho, vY)};')
        lines.append(f'  log_k_f += log(fmax({vsmall}, cTBC));')
        lines.append(f'  fpdtype_t log_rp = log_k_f + {self.fwd_rate_expr(vlog_cs)};')
        lines.append('  fpdtype_t rp = exp(log_rp);')
        if self.reversible:
            lines.append(f'  {self.reverse_rate_block(vlog_cs, vgbs, vprefRuT, vprefRuTinv)}')
        lines.append(f'  {self.accumulate_omega_expr(vomega)}')
        lines.append('}')
        return '\n'.join(lines)


    @clean_csigns
    def jac_block(self, nvars, vsmall='1e-300'):
        lines = [f'// R{self.index}: {self.equation} (Jacobian)']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr()};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr("rho", "q")};')
        lines.append('  log_k_f += log(cTBC);')
        lines += [f'  {l}' for l in self._jac_rates_lines(need_rp=True)]
        lines += [f'  {l}' for l in self._jac_dT_lines(nvars)]
        lines += [f'  {l}' for l in self._jac_dC_lines(nvars, vsmall,
                                                       tb_factor='1.0/cTBC')]
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
                   vprefRuT='log_prefRuT', vprefRuTinv='log_prefRuTinv',
                   vsmall='1e-300'):
        lines = [f'// R{self.index}: {self.equation} (Lindemann)']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr(vlogT, vTinv)};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr(vrho, vY)};')
        lines.append(f'  fpdtype_t log_Pr = log(fmax({vsmall}, cTBC)) + {self.log_pr_expr(vlogT, vTinv)};')
        lines.append('  fpdtype_t log_pmod = -log1p(exp(-log_Pr));')
        lines.append('  log_k_f += log_pmod;')
        lines.append(f'  fpdtype_t log_rp = log_k_f + {self.fwd_rate_expr(vlog_cs)};')
        lines.append('  fpdtype_t rp = exp(log_rp);')
        if self.reversible:
            lines.append(f'  {self.reverse_rate_block(vlog_cs, vgbs, vprefRuT, vprefRuTinv)}')
        lines.append(f'  {self.accumulate_omega_expr(vomega)}')
        lines.append('}')
        return '\n'.join(lines)


    def dlog_pr_dT_expr(self, vTinv='Tinv'):
        r, lo = self.rate, self.low_rate
        return f'({lo.b - r.b} + {lo.Ea - r.Ea}*{vTinv})*{vTinv}'

    @clean_csigns
    def jac_block(self, nvars, vsmall='1e-300'):
        lines = [f'// R{self.index}: {self.equation} (Lindemann, Jacobian)']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr()};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr("rho", "q")};')
        lines.append(f'  fpdtype_t log_Pr = log(cTBC) + {self.log_pr_expr()};')
        lines.append('  fpdtype_t log_pmod = -log1p(exp(-log_Pr));')
        lines.append('  log_k_f += log_pmod;')
        lines.append('  fpdtype_t dpm = 1.0/(1.0 + exp(log_Pr));')
        lines.append(f'  fpdtype_t dlPr_dT = {self.dlog_pr_dT_expr()};')
        lines += [f'  {l}' for l in self._jac_rates_lines(need_rp=True)]
        lines += [f'  {l}' for l in self._jac_dT_lines(nvars,
                                                       extra_dT='dpm*dlPr_dT')]
        lines += [f'  {l}' for l in self._jac_dC_lines(nvars, vsmall,
                                                       tb_factor='dpm/cTBC')]
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

    def _dlog10fcent_dT_block(self, vsmall, vT='T', vTinv='Tinv'):
        # dL/dT with L = log10(Fcent), via direct evaluation of Fcent and
        # its temperature derivative
        alpha, Tsss, Ts = self.troe_A, self.troe_T3, self.troe_T1
        Tss = self.troe_T2
        lines = ['{ fpdtype_t _Fc = 0.0, _dFc = 0.0;']
        if alpha != 1.0:
            lines.append(f'  {{ fpdtype_t _e = exp(-{vT}*{1.0/Tsss}); '
                         f'_Fc += {1.0 - alpha}*_e; '
                         f'_dFc += {-(1.0 - alpha)/Tsss}*_e; }}')
        if alpha != 0.0:
            lines.append(f'  {{ fpdtype_t _e = exp(-{vT}*{1.0/Ts}); '
                         f'_Fc += {alpha}*_e; '
                         f'_dFc += {-alpha/Ts}*_e; }}')
        if Tss is not None and Tss != 0.0:
            lines.append(f'  {{ fpdtype_t _e = exp(-{Tss}*{vTinv}); '
                         f'_Fc += _e; '
                         f'_dFc += {Tss}*{vTinv}*{vTinv}*_e; }}')
        lines.append(f'  dL_dT = _dFc/(fmax(_Fc, {vsmall})*{math.log(10.0)});')
        lines.append('}')
        return '\n  '.join(lines)

    @clean_csigns
    def jac_block(self, nvars, vsmall='1e-300'):
        ln10 = math.log(10.0)
        da_dL = -0.67
        dd_dL = -1.27 + 0.14*0.67
        lines = [f'// R{self.index}: {self.equation} (Troe, Jacobian)']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr()};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr("rho", "q")};')
        lines.append(f'  fpdtype_t log_Pr = log(cTBC) + {self.log_pr_expr()};')
        lines.append(f'  {self._log10_fcent_block("T", "Tinv")}')
        lines.append(f'  fpdtype_t C = -0.4 - 0.67*log10Fcent;')
        lines.append(f'  fpdtype_t N = 0.75 - 1.27*log10Fcent;')
        lines.append(f'  fpdtype_t log10_Pr = log_Pr * {1.0/ln10};')
        lines.append(f'  fpdtype_t A_troe = log10_Pr + C;')
        lines.append(f'  fpdtype_t D_troe = N - 0.14*A_troe;')
        lines.append('  fpdtype_t _s = D_troe*D_troe + A_troe*A_troe;')
        lines.append('  fpdtype_t _sinv = 1.0/_s;')
        lines.append('  fpdtype_t _g = D_troe*D_troe*_sinv;')
        lines.append(f'  fpdtype_t log_F_pdr = log10Fcent*_g * {ln10};')
        lines.append('  fpdtype_t log_pmod = -log1p(exp(-log_Pr)) + log_F_pdr;')
        lines.append('  log_k_f += log_pmod;')
        lines.append('  fpdtype_t dpm = 1.0/(1.0 + exp(log_Pr));')
        lines.append(f'  fpdtype_t dlPr_dT = {self.dlog_pr_dT_expr()};')
        # f1*g^2*df1_dA = A*D*N/(A^2 + D^2)^2 etc.: grouped so every
        # intermediate stays bounded as D_troe -> 0
        lines.append('  fpdtype_t _ad = A_troe*D_troe*_sinv*_sinv;')
        lines.append(f'  fpdtype_t _dlF_dA = -{2.0*ln10}*log10Fcent'
                     f'*_ad*N;')
        lines.append(f'  fpdtype_t _dlF_dL = {ln10}*(_g '
                     f'- 2.0*log10Fcent*_ad*(({da_dL})*D_troe '
                     f'- A_troe*({dd_dL})));')
        lines.append('  fpdtype_t dL_dT;')
        lines.append(f'  {self._dlog10fcent_dT_block(vsmall)}')
        lines += [f'  {l}' for l in self._jac_rates_lines(need_rp=True)]
        extra_dT = (f'dpm*dlPr_dT + _dlF_dA*dlPr_dT*{1.0/ln10} '
                    f'+ _dlF_dL*dL_dT')
        lines += [f'  {l}' for l in self._jac_dT_lines(nvars,
                                                       extra_dT=extra_dT)]
        tb_factor = f'(dpm + _dlF_dA*{1.0/ln10})/cTBC'
        lines += [f'  {l}' for l in self._jac_dC_lines(nvars, vsmall,
                                                       tb_factor=tb_factor)]
        lines.append('}')
        return '\n'.join(lines)

    @clean_csigns
    def rate_block(self, vlogT='logT', vTinv='Tinv', vlog_cs='log_cs',
                   vgbs='gbs', vrho='rho', vY='Y', vomega='omega',
                   vprefRuT='log_prefRuT', vprefRuTinv='log_prefRuTinv',
                   vT='T', vsmall='1e-300'):
        ln10 = math.log(10.0)
        lines = [f'// R{self.index}: {self.equation} (Troe)']
        lines.append('{')
        lines.append(f'  fpdtype_t log_k_f = {self.log_kf_expr(vlogT, vTinv)};')
        lines.append(f'  fpdtype_t cTBC = {self.ctbc_expr(vrho, vY)};')
        lines.append(f'  fpdtype_t log_Pr = log(fmax({vsmall}, cTBC)) + {self.log_pr_expr(vlogT, vTinv)};')
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
