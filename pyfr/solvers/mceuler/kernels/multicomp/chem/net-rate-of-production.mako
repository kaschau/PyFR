<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%! import math %>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

## Log-space Arrhenius: ln k = ln A + b ln T - (Ea/Ru)/T
<%def name="arrh(r)" filter="trim">\
${r.log_A}${f' + ({r.b}*logT)' if r.b else ''}${f' - ({r.Ea}*Tinv)' if r.Ea else ''}\
</%def>

## Forward concentration sum: sum_n nu_n log c_n (or explicit orders)
<%def name="fwd_sum(rxn)" filter="trim">\
${' + '.join(f'({v})*log_cs[{n}]' for n, v in rxn.fwd_exps) or '0.0'}\
</%def>

## Third-body concentration: default*ctot + sum_k (eff_k - default)*c_k
<%def name="ctbc(rxn)" filter="trim">\
<%
    parts = []
    if rxn.default_efficiency == 1.0:
        parts.append('ctot')
    elif rxn.default_efficiency != 0.0:
        parts.append(f'{rxn.default_efficiency}*ctot')
    if rxn.tbc_devs:
        parts.append(' + '.join(f'({v})*u[{n}]' for n, v in rxn.tbc_devs))
%>\
${' + '.join(parts) or '0.0'}\
</%def>

## Exponent ceiling: the largest rate of progress whose omega
## contributions provably stay finite, rp <= fpdtype_max
## / (2*max_n MW_n*sum_i |nu_in|); the factor 2 covers accumulation
## rounding
<%def name="expclamp()" filter="trim">\
${math.log(fpdtype_max) - math.log(2.0*mcf.max_omega_gain)}\
</%def>

## Net rate of progress from the affinity ln(rop_f/rop_r) =
## -sum_n nu_n (gbs_n + ln c_n), less the standard-state pressure term.
## Since that affinity is by definition log_rp - log_rp_r, the reverse
## exponent is never accumulated separately: the surviving direction is
## log_rp - fmin(rp_diff, 0), and expm1 preserves the net rate near
## equilibrium, where the two directions cancel
<%def name="rev_net_rate(rxn)" filter="trim">\
<%
    parts = []
    # nu_total * ln(p0/(Ru T)); the reverse sign folds in since
    # ln(Ru T/p0) = -ln(p0/(Ru T)), so both directions use log_prefRuT
    if rxn.nu_total != 0:
        parts.append(f'{rxn.nu_total}*log_prefRuT')
    parts += [f'- ({v})*(gbs[{n}] + log_cs[{n}])' for n, v in rxn.nu_terms]
%>\
    fpdtype_t rp_diff = ${' '.join(parts)};
    fpdtype_t rp = copysign(exp(fmin(log_rp - fmin(rp_diff, 0.0),
                                     ${expclamp()})),
                            rp_diff)*(-expm1(-fabs(rp_diff)));\
</%def>

<%def name="fwd_only_rate()" filter="trim">\
    fpdtype_t rp = exp(fmin(log_rp, ${expclamp()}));\
</%def>

<%def name="omega_accum(rxn)" filter="trim">\
% for n, v in rxn.nu_terms:
    omega[${n}] += ${v*mcf[n].MW}*rp;
% endfor
</%def>

## log(Pr/(1 + Pr)) in softplus form: exact and overflow-proof for any
## log_Pr, unlike -log1p(exp(-log_Pr))
<%def name="log_pmod_expr()" filter="trim">\
fmin(log_Pr, 0.0) - log1p(exp(-fabs(log_Pr)))\
</%def>

## log10(Fcent) via signed log-sum-exp about the largest term; when
## negative terms can take Fcent through zero the difference is clamped
## like Cantera's log10(max(Fcent, small))
<%def name="troe_fcent(rxn)" filter="trim">\
<%
    def texpr(t):
        pos, c0, cT, cTinv = t
        e = ([f'{c0}'] if c0 else []) \
          + ([f'({cT})*T'] if cT else []) \
          + ([f'({cTinv})*Tinv'] if cTinv else [])
        return ' + '.join(e) or '0.0'

    terms = rxn.fcent_terms
    pos = [texpr(t) for t in terms if t[0]]
    neg = [texpr(t) for t in terms if not t[0]]
    log10e = math.log10(math.e)
%>\
% if not pos:
    fpdtype_t log10Fcent = ${math.log10(fpdtype_min)};
% elif len(terms) == 1:
    fpdtype_t log10Fcent = (${pos[0]}) * ${log10e};
% else:
<%
    refex = texpr(terms[0])
    for t in terms[1:]:
        refex = f'fmax({texpr(t)}, {refex})'
    s = ' + '.join(f'exp({e} - fref)' for e in pos)
    if neg:
        s += ' - ' + ' - '.join(f'exp({e} - fref)' for e in neg)
        s = f'fmax({fpdtype_min}, {s})'
%>\
    fpdtype_t log10Fcent;
    { fpdtype_t fref = ${refex};
      log10Fcent = (fref + log(${s})) * ${log10e}; }
% endif
</%def>

<%def name="elementary_block(rxn)" filter="trim">\
  // R${rxn.index}: ${rxn.equation}
  {
    fpdtype_t log_k_f = ${arrh(rxn.rate)};
    fpdtype_t log_rp = log_k_f + ${fwd_sum(rxn)};
${rev_net_rate(rxn) if rxn.reversible else fwd_only_rate()}
${omega_accum(rxn)}
  }
</%def>

<%def name="three_body_block(rxn)" filter="trim">\
  // R${rxn.index}: ${rxn.equation}
  {
    fpdtype_t log_k_f = ${arrh(rxn.rate)};
    fpdtype_t cTBC = ${ctbc(rxn)};
    log_k_f += log(fmax(${fpdtype_min}, cTBC));
    fpdtype_t log_rp = log_k_f + ${fwd_sum(rxn)};
${rev_net_rate(rxn) if rxn.reversible else fwd_only_rate()}
${omega_accum(rxn)}
  }
</%def>

<%def name="lindemann_block(rxn)" filter="trim">\
  // R${rxn.index}: ${rxn.equation} (Lindemann)
  {
    fpdtype_t log_k_f = ${arrh(rxn.rate)};
    fpdtype_t cTBC = ${ctbc(rxn)};
    fpdtype_t log_Pr = log(fmax(${fpdtype_min}, cTBC)) + ${arrh(rxn.pr_rate)};
    log_k_f += ${log_pmod_expr()};
    fpdtype_t log_rp = log_k_f + ${fwd_sum(rxn)};
${rev_net_rate(rxn) if rxn.reversible else fwd_only_rate()}
${omega_accum(rxn)}
  }
</%def>

<%def name="troe_block(rxn)" filter="trim">\
  // R${rxn.index}: ${rxn.equation} (Troe)
  {
    fpdtype_t log_k_f = ${arrh(rxn.rate)};
    fpdtype_t cTBC = ${ctbc(rxn)};
    fpdtype_t log_Pr = log(fmax(${fpdtype_min}, cTBC)) + ${arrh(rxn.pr_rate)};
${troe_fcent(rxn)}
    fpdtype_t C = -0.4 - 0.67*log10Fcent;
    fpdtype_t N = 0.75 - 1.27*log10Fcent;
    fpdtype_t A_troe = log_Pr*${1.0/math.log(10.0)} + C;
    fpdtype_t f1 = A_troe/(N - 0.14*A_troe);
    fpdtype_t log_F_pdr = log10Fcent/(1.0 + f1*f1) * ${math.log(10.0)};
    log_k_f += ${log_pmod_expr()} + log_F_pdr;
    fpdtype_t log_rp = log_k_f + ${fwd_sum(rxn)};
${rev_net_rate(rxn) if rxn.reversible else fwd_only_rate()}
${omega_accum(rxn)}
  }
</%def>

<%pyfr:macro name='net_rate_of_production' params='u, T, omega'>
## Log-concentration floor for absent species, derived from the
## mechanism's own rate constants and affinities over its thermo fit
## range (see Chemistry.log_c_floor). Kept as shallow as the mechanism
## allows: low-precision builds round the rate exponent at one ulp of
## its largest intermediate. fmax also absorbs the -inf/NaN from log
## of zero/negative concentrations
<% log_cfloor = mcf.log_c_floor %>

  % for n in range(mcf.ns):
    omega[${n}] = 0.0;
  % endfor

  // Concentrations (log space)
  fpdtype_t log_cs[${mcf.ns}];
  % for n in mcf.conc_species:
    log_cs[${n}] = fmax(${log_cfloor}, log(u[${n}]*${1.0/mcf[n].MW}));
  % endfor

  fpdtype_t logT = log(T);
  fpdtype_t Tinv = 1.0/T;
% if mcf.has_reversible:
  // Gibbs free energy (non-dimensional G/RT)
  fpdtype_t gbs[${mcf.ns}];
  fpdtype_t log_prefRuT = ${math.log(101325.0/mcf.Ru)} - logT;
  % for n in mcf.gbs_species:
    gbs[${n}] = ${mcf[n].gbs_expr('T', 'logT', 'Tinv')};
  % endfor
% endif
% if mcf.has_third_body:
  // Total molar concentration for third-body efficiencies
  fpdtype_t ctot = ${' + '.join(f'({c})*u[{n}]' for n, c in mcf.ctot_coeffs)};
% endif

  // Elementary reactions
  % for rxn in mcf.reactions_by_type('elementary'):
    % if not rxn.zero_rate:
${elementary_block(rxn)}
    % endif
  % endfor

  // Three-body reactions
  % for rxn in mcf.reactions_by_type('three-body'):
    % if not rxn.zero_rate:
${three_body_block(rxn)}
    % endif
  % endfor

  // Lindemann falloff reactions
  % for rxn in mcf.reactions_by_type('falloff-lindemann'):
    % if not rxn.zero_rate:
${lindemann_block(rxn)}
    % endif
  % endfor

  // Troe falloff reactions
  % for rxn in mcf.reactions_by_type('falloff-troe'):
    % if not rxn.zero_rate:
${troe_block(rxn)}
    % endif
  % endfor

</%pyfr:macro>
