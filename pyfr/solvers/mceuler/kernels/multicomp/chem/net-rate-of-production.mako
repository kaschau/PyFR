<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.multicomp.makoutil' name='mc'/>
<% import math %>\

<% ns, vix, Eix, rhoix, pix, Tix = mc.thermix(c['ns'], ndims) %>\

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

<%def name="logSumExp2(log_a, log_b, is_diff)">
  {
    // Two-term log-sum-exp or log-diff-exp
    fpdtype_t log_ref = fmax(${log_a}, ${log_b});
  % if is_diff:
    log_result = log_ref + log(exp(${log_a} - log_ref) - exp(${log_b} - log_ref));
  % else:
    log_result = log_ref + log(exp(${log_a} - log_ref) + exp(${log_b} - log_ref));
  % endif
  }
</%def>\

<%def name="logSumExp3(log_a, log_b, log_c)">
  {
    // Three-term log-sum-exp
    fpdtype_t log_ref = fmax(fmax(${log_a}, ${log_b}), ${log_c});
    log_result = log_ref + log(exp(${log_a} - log_ref) + exp(${log_b} - log_ref) + exp(${log_c} - log_ref));
  }
</%def>\

<%def name="computeLog10Fcent(alpha, Tsss, Ts, Tss=0.0)">
  <%
  alpha = float(alpha)
  Tsss = float(Tsss)
  Ts = float(Ts)
  Tss = float(Tss)
  ln_to_log10 = 1.0 / math.log(10.0)  # Convert natural log to log10
  is_three_param = (Tss == 0.0)
  %>
  % if alpha == 0.0:
    % if is_three_param:
      // α = 0, 3-param: F_cent = exp(-T/T***)
      fpdtype_t log10Fcent = -T*${ln_to_log10/Tsss};
    % else:
      // α = 0, 4-param: F_cent = exp(-T/T***) + exp(-T**/T)
      fpdtype_t log_result;
      ${logSumExp2(f"-T*{1.0/Tsss}", f"-{Tss}*Tinv", False)}
      fpdtype_t log10Fcent = log_result * ${ln_to_log10};
    % endif
  % elif alpha == 1.0:
    % if is_three_param:
      // α = 1, 3-param: F_cent = exp(-T/T*)
      fpdtype_t log10Fcent = -T*${ln_to_log10/Ts};
    % else:
      // α = 1, 4-param: F_cent = exp(-T/T*) + exp(-T**/T)
      fpdtype_t log_result;
      ${logSumExp2(f"-T*{1.0/Ts}", f"-{Tss}*Tinv", False)}
      fpdtype_t log10Fcent = log_result * ${ln_to_log10};
    % endif
  % elif alpha > 1.0:
    % if is_three_param:
      // α > 1, 3-param: F_cent = α*exp(-T/T*) - (α-1)*exp(-T/T***)
      fpdtype_t log_result;
      ${logSumExp2(f"{math.log(alpha)} - T*{1.0/Ts}", f"{math.log(alpha - 1.0)} - T*{1.0/Tsss}", True)}
      fpdtype_t log10Fcent = log_result * ${ln_to_log10};
    % else:
      // α > 1, 4-param: F_cent = α*exp(-T/T*) - (α-1)*exp(-T/T***) + exp(-T**/T)
      // Compute (α*exp(-T/T*) + exp(-T**/T)) - (α-1)*exp(-T/T***)
      fpdtype_t log_result;
      ${logSumExp2(f"{math.log(alpha)} - T*{1.0/Ts}", f"-{Tss}*Tinv", False)}
      fpdtype_t log_pos_sum = log_result;
      ${logSumExp2("log_pos_sum", f"{math.log(alpha - 1.0)} - T*{1.0/Tsss}", True)}
      fpdtype_t log10Fcent = log_result * ${ln_to_log10};
    % endif
  % elif alpha > 0.0:
    % if is_three_param:
      // 0 < α < 1, 3-param: F_cent = (1-α)*exp(-T/T***) + α*exp(-T/T*)
      fpdtype_t log_result;
      ${logSumExp2(f"{math.log(1.0 - alpha)} - T*{1.0/Tsss}", f"{math.log(alpha)} - T*{1.0/Ts}", False)}
      fpdtype_t log10Fcent = log_result * ${ln_to_log10};
    % else:
      // 0 < α < 1, 4-param: F_cent = (1-α)*exp(-T/T***) + α*exp(-T/T*) + exp(-T**/T)
      fpdtype_t log_result;
      ${logSumExp3(f"{math.log(1.0 - alpha)} - T*{1.0/Tsss}", f"{math.log(alpha)} - T*{1.0/Ts}", f"-{Tss}*Tinv")}
      fpdtype_t log10Fcent = log_result * ${ln_to_log10};
    % endif
  % else:
    % if is_three_param:
      // α < 0, 3-param: F_cent = (1-α)*exp(-T/T***) - |α|*exp(-T/T*)
      fpdtype_t log_result;
      ${logSumExp2(f"{math.log(1.0 - alpha)} - T*{1.0/Tsss}", f"{math.log(-alpha)} - T*{1.0/Ts}", True)}
      fpdtype_t log10Fcent = log_result * ${ln_to_log10};
    % else:
      // α < 0, 4-param: F_cent = (1-α)*exp(-T/T***) - |α|*exp(-T/T*) + exp(-T**/T)
      // Compute ((1-α)*exp(-T/T***) + exp(-T**/T)) - |α|*exp(-T/T*)
      fpdtype_t log_result;
      ${logSumExp2(f"{math.log(1.0 - alpha)} - T*{1.0/Tsss}", f"-{Tss}*Tinv", False)}
      fpdtype_t log_pos_sum = log_result;
      ${logSumExp2("log_pos_sum", f"{math.log(-alpha)} - T*{1.0/Ts}", True)}
      fpdtype_t log10Fcent = log_result * ${ln_to_log10};
    % endif
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
      gbs[${n}] = ${mc.nasa_gbs(fast_coeff[n])};
    % else:
      if (T < ${T_cutoff[n]}){
        gbs[${n}] = ${mc.nasa_gbs(NASA7_Tlow[n])};
      }else{
        gbs[${n}] = ${mc.nasa_gbs(NASA7_Thigh[n])};
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
    ${computeLog10Fcent(alpha, Tsss, Ts, Tss)}
    fpdtype_t C = -0.4 - 0.67*log10Fcent;
    fpdtype_t N = 0.75 - 1.27*log10Fcent;
    fpdtype_t log_Pr = log_cTBC + ${logRateConst(A_o[i]/A_f[i], m_o[i]-m_f[i], Ea_o[i]-Ea_f[i])};
    fpdtype_t log10_Pr = log_Pr * ${1.0 / math.log(10.0)}; // Convert ln to log10
    fpdtype_t A = log10_Pr + C;
    fpdtype_t f1 = A/(N - 0.14*A);
    fpdtype_t log_F_pdr = log10Fcent/(1.0+f1*f1) * ${math.log(10.0)}; // Convert log10 to ln
    fpdtype_t log_pmod = log_Pr - log(1.0 + exp(log_Pr)) + log_F_pdr;
    log_k_f += log_pmod;
  % elif c['r_type'][i] == 'SRI':
  <% raise NotImplementedError("SRI reactions not supporeted")%>
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