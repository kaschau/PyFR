<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% import math %>\

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\

<% nr = c['Ea_f'].shape[0] %>\
<% MW = c['MW'] %>\
<% Ru = c['Ru'] %>\
<% fast_props = 'fast_coeff' in c %>\
% if fast_props:
<% fast_coeff = c['fast_coeff'] %>\
% else:
<% T_cutoff = c['T_cutoff'] %>\
<% NASA7_Thigh = c['NASA7_Thigh'] %>\
<% NASA7_Tlow = c['NASA7_Tlow'] %>\
% endif\

<%def name="logRateConst(A, m, Ea)">
  ${math.log(A)}+(${m}*logT)-(${Ea}*Tinv)
</%def>\

<%def name="Kc_log(nusum)">
  <%
  nusum = float(nusum)
  if nusum > 0.0:
      log_term = f"{nusum}*log_prefRuTinv"
  elif nusum < 0.0:
      log_term = f"{-nusum}*log_prefRuT"
  else:
      log_term = "0.0"
  %>
  -(log_Kp + ${log_term})
</%def>\

<%def name="troeThreeParam(alpha, Tsss, Ts)">
  <%
  alpha = float(alpha)
  Tsss = float(Tsss)
  Ts = float(Ts)
  %>
  % if alpha == 0.0:
  // Special case: α = 0, F_cent = exp(-T/T***)
  fpdtype_t log10Fcent = -T*${0.4342944819032518/Tsss}; // ln to log10
  % elif alpha == 1.0:
  // Special case: α = 1, F_cent = exp(-T/T*)
  fpdtype_t log10Fcent = -T*${0.4342944819032518/Ts}; // ln to log10
  % elif alpha > 1.0:
  // Alpha > 1: F_cent = α*exp(-T/T*) - (α-1)*exp(-T/T***)
  fpdtype_t log_term1 = ${math.log(alpha)} - T*${1.0/Ts};
  fpdtype_t log_term2 = ${math.log(alpha - 1.0)} - T*${1.0/Tsss};
  // Use log-difference with fmax reference selection
  fpdtype_t log_ref = fmax(log_term1, log_term2);
  fpdtype_t log_sum = log_ref + log(exp(log_term1 - log_ref) - exp(log_term2 - log_ref));
  fpdtype_t log10Fcent = log_sum * 0.4342944819032518; // ln to log10
  % elif alpha > 0.0:
  // 0 < alpha < 1: F_cent = (1-α)*exp(-T/T***) + α*exp(-T/T*)
  fpdtype_t log_term1 = ${math.log(1.0 - alpha)} - T*${1.0/Tsss};
  fpdtype_t log_term2 = ${math.log(alpha)} - T*${1.0/Ts};
  // Use log-sum-exp with fmax reference selection
  fpdtype_t log_ref = fmax(log_term1, log_term2);
  fpdtype_t log_sum = log_ref + log(exp(log_term1 - log_ref) + exp(log_term2 - log_ref));
  fpdtype_t log10Fcent = log_sum * 0.4342944819032518; // ln to log10
  % else:
  // Negative alpha: F_cent = (1-α)*exp(-T/T***) - |α|*exp(-T/T*)
  fpdtype_t log_term1 = ${math.log(1.0 - alpha)} - T*${1.0/Tsss};
  fpdtype_t log_term2 = ${math.log(-alpha)} - T*${1.0/Ts};
  // Use log-difference: log(a - b) = log_ref + log(exp(log_a - log_ref) - exp(log_b - log_ref))
  fpdtype_t log_ref = fmax(log_term1, log_term2);
  fpdtype_t log_sum = log_ref + log(exp(log_term1 - log_ref) - exp(log_term2 - log_ref));
  fpdtype_t log10Fcent = log_sum * 0.4342944819032518; // ln to log10
  % endif
</%def>\

<%def name="troeFourParam(alpha, Tsss, Ts, Tss)">
  <%
  alpha = float(alpha)
  Tsss = float(Tsss)
  Ts = float(Ts)
  Tss = float(Tss)
  %>
  % if alpha == 0.0:
  // Special case: α = 0, F_cent = exp(-T/T***) + exp(-T**/T)
  fpdtype_t log_term1 = -T*${1.0/Tsss};
  fpdtype_t log_term3 = -${Tss}*Tinv;
  fpdtype_t log_ref = fmax(log_term1, log_term3);
  fpdtype_t log_sum = log_ref + log(exp(log_term1 - log_ref) + exp(log_term3 - log_ref));
  fpdtype_t log10Fcent = log_sum * 0.4342944819032518; // ln to log10
  % elif alpha == 1.0:
  // Special case: α = 1, F_cent = exp(-T/T*) + exp(-T**/T)
  fpdtype_t log_term2 = -T*${1.0/Ts};
  fpdtype_t log_term3 = -${Tss}*Tinv;
  fpdtype_t log_ref = fmax(log_term2, log_term3);
  fpdtype_t log_sum = log_ref + log(exp(log_term2 - log_ref) + exp(log_term3 - log_ref));
  fpdtype_t log10Fcent = log_sum * 0.4342944819032518; // ln to log10
  % elif alpha > 1.0:
  // Alpha > 1: F_cent = α*exp(-T/T*) - |1-α|*exp(-T/T***) + exp(-T**/T)
  fpdtype_t log_term1 = ${math.log(alpha)} - T*${1.0/Ts};
  fpdtype_t log_term2 = ${math.log(alpha - 1.0)} - T*${1.0/Tsss};
  fpdtype_t log_term3 = -${Tss}*Tinv;
  // F_cent = term1 + term3 - term2 = (term1 + term3) - term2
  // First compute log(term1 + term3) using log-sum-exp
  fpdtype_t log_pos_ref = fmax(log_term1, log_term3);
  fpdtype_t log_pos_sum = log_pos_ref + log(exp(log_term1 - log_pos_ref) + exp(log_term3 - log_pos_ref));
  // Then compute log(pos_sum - term2) using log-difference
  fpdtype_t log_ref = fmax(log_pos_sum, log_term2);
  fpdtype_t log_sum = log_ref + log(exp(log_pos_sum - log_ref) - exp(log_term2 - log_ref));
  fpdtype_t log10Fcent = log_sum * 0.4342944819032518; // ln to log10
  % elif alpha > 0.0:
  // 0 < alpha < 1: F_cent = (1-α)*exp(-T/T***) + α*exp(-T/T*) + exp(-T**/T)
  // Choose reference to avoid single precision overflow (exp arg > 88)
  fpdtype_t log_term1 = ${math.log(1.0 - alpha)} - T*${1.0/Tsss};
  fpdtype_t log_term2 = ${math.log(alpha)} - T*${1.0/Ts};
  fpdtype_t log_term3 = -${Tss}*Tinv;
  fpdtype_t log_ref = fmax(fmax(log_term1, log_term2), log_term3);
  fpdtype_t log_sum = log_ref + log(exp(log_term1 - log_ref) + exp(log_term2 - log_ref) + exp(log_term3 - log_ref));
  fpdtype_t log10Fcent = log_sum * 0.4342944819032518; // ln to log10
  % else:
  // Negative alpha: F_cent = (1-α)*exp(-T/T***) - |α|*exp(-T/T*) + exp(-T**/T)
  fpdtype_t log_term1 = ${math.log(1.0 - alpha)} - T*${1.0/Tsss};
  fpdtype_t log_term2 = ${math.log(-alpha)} - T*${1.0/Ts};  // |α| = -α since α < 0
  fpdtype_t log_term3 = -${Tss}*Tinv;
  // F_cent = term1 + term3 - term2 = (term1 + term3) - term2
  // First compute log(term1 + term3) using log-sum-exp
  fpdtype_t log_pos_ref = fmax(log_term1, log_term3);
  fpdtype_t log_pos_sum = log_pos_ref + log(exp(log_term1 - log_pos_ref) + exp(log_term3 - log_pos_ref));
  // Then compute log(pos_sum - term2) using log-difference
  fpdtype_t log_ref = fmax(log_pos_sum, log_term2);
  fpdtype_t log_sum = log_ref + log(exp(log_pos_sum - log_ref) - exp(log_term2 - log_ref));
  fpdtype_t log10Fcent = log_sum * 0.4342944819032518; // ln to log10
  % endif
</%def>\

<%pyfr:macro name='net_rate_of_production' params='Y, T, rho, omega'>

  % for n in range(ns):
    omega[${n}] = 0.0;  // omega must start at zero
  % endfor
  // Concentrations (log space only)
  fpdtype_t log_cs[${ns}];
  % for n in range(ns):
    log_cs[${n}] = log(fmax(0.0, rho*Y[${n}]*${1.0/c['MW'][n]}));
  % endfor

  // Gibbs energy (kept in log space)
  fpdtype_t gbs[${ns}];
  fpdtype_t logT = log(T);
  fpdtype_t Tinv = 1.0/T;
  fpdtype_t log_prefRuT = ${math.log(101325.0/c['Ru'])} + log(Tinv);
  fpdtype_t log_prefRuTinv = ${math.log(c['Ru']/101325.0)} + logT;
  % for n in range(ns):
    // ${c['names'][n]} Properties
    % if fast_props:
      gbs[${n}] = ${pyfr.nasa_gbs(fast_coeff[n])};
    % else:
      if (T < ${T_cutoff[n]}){
        gbs[${n}] = ${pyfr.nasa_gbs(NASA7_Tlow[n])};
      }else{
        gbs[${n}] = ${pyfr.nasa_gbs(NASA7_Thigh[n])};
      }
    % endif
  % endfor

  // Rate constants, Falloff Mods, new rates of progress
  <% A_f = c['A_f'] %>\
  <% m_f = c['m_f'] %>\
  <% Ea_f = c['Ea_f'] %>\
  <% A_o = c['A_o'] %>\
  <% m_o = c['m_o'] %>\
  <% Ea_o = c['Ea_o'] %>\
  <% nu_f = c['nu_f'] %>\
  <% nu_b = c['nu_b'] %>\

  % for i in range(nr):
  // Reaction ${i} - ${c['r_type'][i]}
  {
  fpdtype_t log_k_f = ${logRateConst(A_f[i], m_f[i], Ea_f[i])};
  % if sum(c['aij'][i]) > 0.0:
  // Three body reaction
  fpdtype_t cTBC = rho * (${"+".join([f"({eff/c['MW'][n]})*Y[{n}]" for n,eff in enumerate(c['aij'][i]) if eff != 0.0])});
  fpdtype_t log_cTBC = log(cTBC);
  % endif
  % if c['r_type'][i] == 'three-body-Arrhenius':
    log_k_f += log_cTBC;
  % elif c['r_type'][i] == 'falloff-Lindemann':
    // Lindemann Reaction (log space)
    fpdtype_t log_Pr = log_cTBC + ${logRateConst(A_o[i]/A_f[i], m_o[i]-m_f[i], Ea_o[i]-Ea_f[i])};
    fpdtype_t log_pmod = log_Pr - log(1.0 + exp(log_Pr));
    log_k_f += log_pmod;
  % elif c['r_type'][i] == 'falloff-Troe':
    <% alpha = c['fall_coeffs'][i][0]%>\
    <% Tsss = c['fall_coeffs'][i][1]%>\
    <% Ts = c['fall_coeffs'][i][2]%>\
    <% Tss = c['fall_coeffs'][i][3]%>\
    % if Tss == 0.0: #Three Parameter Troe form
      ${troeThreeParam(alpha, Tsss, Ts)}
    % else: # Four Parameter Troe form
      ${troeFourParam(alpha, Tsss, Ts, Tss)}
    % endif
    fpdtype_t C = -0.4 - 0.67*log10Fcent;
    fpdtype_t N = 0.75 - 1.27*log10Fcent;
    fpdtype_t log_Pr = log_cTBC + ${logRateConst(A_o[i]/A_f[i], m_o[i]-m_f[i], Ea_o[i]-Ea_f[i])};
    fpdtype_t log10_Pr = log_Pr * 0.4342944819032518; // log_Pr / ln(10)
    fpdtype_t A = log10_Pr + C;
    fpdtype_t f1 = A/(N - 0.14*A);
    fpdtype_t log_F_pdr = log10Fcent/(1.0+f1*f1) * 2.302585092994046; // ln(10)
    fpdtype_t log_pmod = log_Pr - log(1.0 + exp(log_Pr)) + log_F_pdr;
    log_k_f += log_pmod;
  % elif c['r_type'][i] == 'SRI':
  <% raise ImplementedError("SRI reactions not supporeted")%>
  % endif

  // Set rates of progress (log space)
  <% nu_sum = nu_b[:,i] - nu_f[:,i] %>\
  % if c['r_type'][i] == "Arrhenius Custom Order":
    fpdtype_t log_rp = log_k_f + ${"+".join([f"({v})*log_cs[{n}]" for n,v in enumerate(c['orders'][i]) if float(v) != 0.0])};
  % else:
    fpdtype_t log_rp = log_k_f + ${"+".join([f"({v})*log_cs[{n}]" for n,v in enumerate(nu_f[:,i]) if float(v) != 0.0])};
  % endif
  fpdtype_t rp = exp(log_rp);

  % if c['reversible'][i]:
    // Equilibrium constant
    fpdtype_t log_Kp = ${"+".join([f"({v})*gbs[{n}]" for n,v in enumerate(nu_sum) if float(v) != 0.0])};

    // Work in log space to avoid overflow
    fpdtype_t log_k_r = log_k_f - ${Kc_log(sum(nu_sum))};
    fpdtype_t log_rp_reverse = log_k_r + ${"+".join([f"({v})*log_cs[{n}]" for n,v in enumerate(nu_b[:,i]) if float(v) != 0.0])};
    fpdtype_t rp_reverse = exp(log_rp_reverse);
    rp -= rp_reverse;
  % endif

  // Add this reaction to the sources that use it
  % for n in range(ns):
    <% nu = nu_b[n,i] - nu_f[n,i] %>\
    % if abs(nu) > 0.0:
      omega[${n}] += ${nu}*rp;
    % endif
  % endfor
  }
% endfor ##// End reaction loop

// Convert to mass
% for n in range(ns):
  omega[${n}] *= ${MW[n]};
% endfor
</%pyfr:macro>