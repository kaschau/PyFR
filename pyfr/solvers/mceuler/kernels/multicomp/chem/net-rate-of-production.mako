<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% import math %>\

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\

<% nr = c['Ea_f'].shape[0] %>\
<% MW = c['MW'] %>\
<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% fast_props = N7.shape[1] == 7 %>\

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

<%pyfr:macro name='net_rate_of_production' params='q, T, rho, omega'>

  % for n in range(ns):
    omega[${n}] = 0.0;  // omega must start at zero
  % endfor
  // Concentrations (log space only)
  fpdtype_t log_cs[${ns}];
  % for n in range(ns):
    log_cs[${n}] = log(fmax(0.0, rho*q[${n}]*${1.0/c['MW'][n]}));
  % endfor

  // Gibbs energy (kept in log space)
  fpdtype_t gbs[${ns}];
  fpdtype_t logT = log(T);
  fpdtype_t Tinv = 1.0/T;
  fpdtype_t prefRuT = ${101325.0/c['Ru']}*Tinv;
  fpdtype_t prefRuTinv = ${c['Ru']/101325.0}*T;
  fpdtype_t log_prefRuT = log(prefRuT);
  fpdtype_t log_prefRuTinv = log(prefRuTinv);
  % for n in range(ns):
    // ${c['names'][n]} Properties
    % if fast_props:
      gbs[${n}] = ${pyfr.nasa_gbs(N7[n,:], 0)};
    % else:
      if (T < ${N7[n,0]}){
        gbs[${n}] = ${pyfr.nasa_gbs(N7[n,:], 8)};
      }else{
        gbs[${n}] = ${pyfr.nasa_gbs(N7[n,:], 1)};
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
  fpdtype_t cTBC = rho * (${"+".join([f"({eff})*q[{n}]*{1.0/c['MW'][n]}" for n,eff in enumerate(c['aij'][i]) if eff != 0.0])});
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
      // Three Troe Reaction
      ## double log10Fcent = log10((${1.0 - alpha})*exp(-T*${1.0/Tsss}) + ${alpha}*exp(-T*${1.0/Ts}));
      ## Convert to nat log and simplify
      fpdtype_t log10Fcent = ${1.0/math.log(10)}*(-T*${1.0/Tsss} + log(${1.0 - alpha} + ${alpha}*exp(T*${(-Tsss+Ts)/(Tsss*Ts)})));
    % else: # Four Parameter Troe form
      // Four Troe Reaction
      ## double log10Fcent = log10((${1.0 - alpha})*exp(-T*${1.0/Tsss}) + ${alpha}*exp(-T*${1.0/Ts}) + exp(-${Tss}*Tinv));
      ## Convert to nat log and simplify
      fpdtype_t log10Fcent = ${1.0/math.log(10)}*(-T*${1.0/Tsss} + log(${1.0 - alpha} + ${alpha}*exp(T*${(-Tsss+Ts)/(Tsss*Ts)}) + exp(${-Tss}*Tinv + T*${1.0/Tsss})));
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