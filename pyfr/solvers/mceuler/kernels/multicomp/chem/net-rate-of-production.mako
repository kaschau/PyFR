<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% import math %>\

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\

<% nr = c['Ea_f'].shape[0] %>\
<% MW = c['MW'] %>\
<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% fast_props = N7.shape[1] == 7 %>\

<%def name="rateConst(A, m, Ea)">
  <% m = float(m) %>\
  % if m == 0.0 and Ea == 0.0:
    ${A}
  % elif m == 0.0 and Ea != 0.0:
    exp(${math.log(A)}-(${Ea}*Tinv))
  % elif m.is_integer() and Ea == 0.0:
  %   if m < 0.0:
    ${A}${"".join("*Tinv" for _ in range(int(abs(m))))}
  %   elif m > 0.0:
    ${A}${"".join("*T" for _ in range(int(m)))}
  %   endif
  % elif m != 0.0 and Ea == 0.0:
    exp(${math.log(A)}+(${m}*logT))
  % elif Ea != 0.0:
    exp(${math.log(A)}+(${m}*logT)-(${Ea}*Tinv))
  % endif
</%def>\

<%def name="Kcinv(nusum)">
  <% nusum = float(nusum) %>\
  % if nusum != 0.0:
    % if nusum == 1.0:
      prefRuTinv*Kp
    % elif nusum == -1.0:
      Kp*prefRuT
    % elif nusum.is_integer():
      % if nusum > 0.0:
        Kp*${pyfr.intpow("prefRuTinv", nusum)}
      % else:
        Kp*${pyfr.intpow("prefRuT", -nusum)}
      % endif:
    % else:
      Kp*pow(prefRuTinv,${nusum})
    %   endif
  % else:
      Kp
  % endif
</%def>\

<%pyfr:macro name='net_rate_of_production' params='q, T, rho, omega'>

  // Concentrations
  fpdtype_t cs[${ns}];
  // Kahan summation error compensation for omega accumulation
  fpdtype_t omega_c[${ns}] = {0};
  % for n in range(ns):
    omega[${n}] = 0.0;  // omega must start at zero
  % endfor
  % for n in range(ns):
    cs[${n}] = fmax(0.0, rho*q[${n}]*${1.0/c['MW'][n]});
  % endfor

  // Gibbs energy (kept in log space)
  fpdtype_t gbs[${ns}];
  fpdtype_t logT = log(T);
  fpdtype_t Tinv = 1.0/T;
  fpdtype_t prefRuT = ${101325.0/c['Ru']}*Tinv;
  fpdtype_t prefRuTinv = ${c['Ru']/101325.0}*T;
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
  fpdtype_t k_f = ${rateConst(A_f[i], m_f[i], Ea_f[i])};
  % if sum(c['aij'][i]) > 0.0:
  // Three body reaction with Kahan summation
  fpdtype_t cTBC = 0.0;
  fpdtype_t cTBC_c = 0.0;  // compensation
  % for n, eff in enumerate(c['aij'][i]):
    % if eff != 0.0:
    {
      fpdtype_t y = ${eff}*cs[${n}] - cTBC_c;
      fpdtype_t t = cTBC + y;
      cTBC_c = (t - cTBC) - y;
      cTBC = t;
    }
    % endif
  % endfor
  % endif
  % if c['r_type'][i] == 'three-body-Arrhenius':
    k_f *= cTBC;
  % elif c['r_type'][i] == 'falloff-Lindemann':
    // Lindemann Reaction
    fpdtype_t Pr = cTBC*${rateConst(A_o[i]/A_f[i], m_o[i]-m_f[i], Ea_o[i]-Ea_f[i])}; // <- ratio k0/k_f
    fpdtype_t pmod = Pr/(1.0 + Pr);
    k_f *= pmod;
  % elif c['r_type'][i] == 'falloff-Troe':
    <% alpha = c['fall_coeffs'][i][0]%>\
    <% Tsss = c['fall_coeffs'][i][1]%>\
    <% Ts = c['fall_coeffs'][i][2]%>\
    <% Tss = c['fall_coeffs'][i][3]%>\
    % if Tss == 0.0: #Three Parameter Troe form
      // Three Troe Reaction
      fpdtype_t log10Fcent = log10((${1.0 - alpha})*exp(-T*${1.0/Tsss}) + ${alpha}*exp(-T*${1.0/Ts}));
    % else: # Four Parameter Troe form
      // Four Troe Reaction
      fpdtype_t log10Fcent = log10((${1.0 - alpha})*exp(-T*${1.0/Tsss}) + ${alpha}*exp(-T*${1.0/Ts}) + exp(-${Tss}*Tinv));
    % endif
    fpdtype_t C = -0.4 - 0.67*log10Fcent;
    fpdtype_t N = 0.75 - 1.27*log10Fcent;
    fpdtype_t Pr = cTBC*${rateConst(A_o[i]/A_f[i], m_o[i]-m_f[i], Ea_o[i]-Ea_f[i])}; // <- ratio k0/k_f
    fpdtype_t A = log10(Pr) + C;
    fpdtype_t f1 = A/(N - 0.14*A);
    fpdtype_t F_pdr = pow(10.0,log10Fcent/(1.0+f1*f1));
    fpdtype_t pmod = Pr/(1.0 + Pr) * F_pdr;
    k_f *= pmod;
  % elif c['r_type'][i] == 'SRI':
  <% raise ImplementedError("SRI reactions not supporeted")%>
  % endif

  // Set rates of progress
  <% nu_sum = nu_b[:,i] - nu_f[:,i] %>\
  % if c['r_type'][i] == "Arrhenius Custom Order":
    fpdtype_t rp = k_f * ${"*".join([pyfr.intpow(f"cs[{n}]",v) for n,v in enumerate(c['orders'][i]) if float(v) != 0.0])};
  % else:
    fpdtype_t rp = k_f * ${"*".join([pyfr.intpow(f"cs[{n}]",v) for n,v in enumerate(nu_f[:,i]) if float(v) != 0.0])};
  % endif

  % if c['reversible'][i]:
    // Equilibrium constant with Kahan summation
    fpdtype_t log_Kp = 0.0;
    fpdtype_t log_Kp_c = 0.0;  // compensation
    % for n, v in enumerate(nu_sum):
      % if float(v) != 0.0:
    {
      fpdtype_t y = ${v}*gbs[${n}] - log_Kp_c;
      fpdtype_t t = log_Kp + y;
      log_Kp_c = (t - log_Kp) - y;
      log_Kp = t;
    }
      % endif
    % endfor
    fpdtype_t Kp = exp(log_Kp);
    fpdtype_t k_r = ${Kcinv(sum(nu_sum))}*k_f;
    rp -= k_r * ${"*".join([pyfr.intpow(f"cs[{n}]",v) for n,v in enumerate(nu_b[:,i]) if float(v) != 0.0])};
  % endif

  // Add this reaction to the sources that use it (Kahan summation)
  % for n in range(ns):
    <% nu = nu_b[n,i] - nu_f[n,i] %>\
    % if abs(nu) > 0.0:
      {
        fpdtype_t y = ${nu}*rp - omega_c[${n}];
        fpdtype_t t = omega[${n}] + y;
        omega_c[${n}] = (t - omega[${n}]) - y;
        omega[${n}] = t;
      }
    % endif
  % endfor
  }
% endfor ##// End reaction loop

// Convert to mass
% for n in range(ns):
  omega[${n}] *= ${MW[n]};
% endfor
</%pyfr:macro>