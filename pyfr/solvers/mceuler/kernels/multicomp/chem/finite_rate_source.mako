<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.mixture_state'/>

<%def name="rateConst(A, m, Ea)">
% if m == 0.0 and Ea == 0.0:
  ${A}
% elif m == 0.0 and Ea != 0.0:
  exp(log(${A})-(${Ea}*Tinv))
% elif abs(int(m) - m) < 1e-14 and Ea == 0.0:
%   if m < 0:
  ${A}${"".join("*Tinv" for _ in range(int(abs(m))))}
%   elif m > 0:
  ${A}${"".join("*T" for _ in range(int(m)))}
%   endif
% elif m != 0.0 and Ea == 0.0:
  exp(log(${A})+(${m}*logT))
% elif Ea !=0:
  exp(log(${A})+(${m}*logT)-(${Ea}*Tinv))
% endif
</%def>

<%def name="eqConst(nusum)">
% if nusum != 0.0:
%   if nusum == 1.0:
  prefRuT*exp(-dG)
%   elif nusum == -1:
  exp(-dG)/prefRuT
%   else:
  pow(prefRuT,${nusum})*exp(-dG)
%   endif
% else:
  exp(-dG)
% endif
</%def>\

<%pyfr:macro name='finite_rate_source' params='t, u, ploc, src'>
<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\
<% nr = c['Ea_f'].shape[0] %>\
<% MW = c['MW'] %>\
<% N7 = c['NASA7'] %>\
<% div = [1.0, 2.0, 3.0, 4.0, 5.0] %>\

  fpdtype_t rho = u[0];

  // Compute thermodynamic properties
  fpdtype_t prims[${nvars+1}];
  fpdtype_t qh[${3+ns}];
  ${pyfr.expand('mixture_state', 'u', 'prims', 'qh')};

  fpdtype_t T = prims[${ndims+1}];

  // Concentrations
  fpdtype_t cs[${ns}];
% for n in range(ns):
  cs[${n}] = rho*prims[${Yix+n}]*${1.0/c['MW'][n]};
% endfor

  // Gibbs energy
  fpdtype_t gbs[${ns}];
  fpdtype_t logT = log(T);
  fpdtype_t Tinv = 1.0/T;
  fpdtype_t prefRuT = ${101325.0/c['Ru']}*Tinv;
  {
% for n in range(ns):
      // ${c['names'][n]} Properties
      if (T < ${N7[n,0]})
      {
<% m = 8 %>\

            fpdtype_t hi = ${f'+ T*('.join(str(c) for c in N7[n,m:m+5]/div)+')'*4} + ${N7[n, m + 5]}*Tinv;
            fpdtype_t scs = ${N7[n, m + 0]} * logT +
                            T*(${f'+ T*('.join(str(c) for c in N7[n,m+1:m+5]/div[0:-1])+')'*3}) + ${N7[n, m + 6]};
            gbs[${n}] = hi - scs;
        }else
        {
<% m = 1 %>\
            fpdtype_t hi = ${f'+ T*('.join(str(c) for c in N7[n,m:m+5]/div)+')'*4} + ${N7[n, m + 5]}*Tinv;
            fpdtype_t scs = ${N7[n, m + 0]} * logT +
                            T*(${f'+ T*('.join(str(c) for c in N7[n,m+1:m+5]/div[0:-1])+')'*3}) + ${N7[n, m + 6]};
            gbs[${n}] = hi - scs;
        }
% endfor
  }

  // Rate constants, Falloff Mods, new rates of progress
<% A_f = c['A_f'] %>\
<% m_f = c['m_f'] %>\
<% Ea_f = c['Ea_f'] %>\
<% A_o = c['A_o'] %>\
<% m_o = c['m_o'] %>\
<% Ea_o = c['Ea_o'] %>\
  fpdtype_t q[${nr}];

% for i in range(nr):
  <% alpha = c['fall_coeffs'][i][0]%>\
  <% Tsss = c['fall_coeffs'][i][1]%>\
  <% Ts = c['fall_coeffs'][i][2]%>\
  <% Tss = c['fall_coeffs'][i][3]%>\
  <% nu_sum = c['nu_b'][:,i] - c['nu_f'][:,i] %>\

  // Reaction ${i} - ${c['r_type'][i]}
  {
  fpdtype_t k_f = ${rateConst(A_f[i], m_f[i], Ea_f[i])};
  fpdtype_t dG = ${"+".join([f"({s}*gbs[{i}])" for i,s in enumerate(nu_sum) if s != 0.0])};
  fpdtype_t K_c = ${eqConst(sum(nu_sum))};
  % if sum(c['aij'][i]) > 0.0:
  // Three body reaction
  fpdtype_t cTBC = ${"+".join([f"({eff}*cs[{j}])" for j,eff in enumerate(c['aij'][i]) if eff != 0.0])};
  % endif
  % if c['r_type'][i] == 'three-body-Arrhenius':
    k_f *= cTBC;
  % elif c['r_type'][i] == 'falloff-Lindemann':
    // Lindemann Reaction
    fpdtype_t Fcent = 1.0;
    fpdtype_t k0 = ${rateConst(A_o[i], m_o[i], Ea_o[i])};
    fpdtype_t Pr = cTBC*k0/k_f;
    fpdtype_t pmod = Pr/(1.0 + Pr);
    k_f *= pmod;
  % elif c['r_type'][i] == 'falloff-Troe':
    // Troe Reaction
    % if Tss == 0: #Three Parameter Troe form
      fpdtype_t Fcent = (1.0 - (${alpha}))*exp(-T/(${Tsss})) + (${alpha})*exp(-T/({Ts}));
    % else: # Four Parameter Troe form
      fpdtype_t Fcent = (1.0 - (${alpha}))*exp(-T/(${Tsss})) + (${alpha})*exp(-T/(${Ts})) + exp(-(${Tss})*Tinv);
    % endif
    fpdtype_t C = -0.4 - 0.67*log10(Fcent);
    fpdtype_t N = 0.75 - 1.27*log10(Fcent);
    fpdtype_t k0 = ${rateConst(A_o[i], m_o[i], Ea_o[i])};
    fpdtype_t Pr = cTBC*k0/k_f;
    fpdtype_t A = log10(Pr) + C;
    fpdtype_t f1 = A/(N - 0.14*A);
    fpdtype_t F_pdr = pow(10.0,log10(Fcent)/(1.0+f1*f1));
    fpdtype_t pmod = Pr/(1.0 + Pr) * F_pdr;
    k_f *= pmod;
  % elif c['r_type'][i] == 'SRI':
  <% raise ImplementedError("SRI reactions not supporeted")%>
  % endif
  fpdtype_t q_f = k_f * ${"*".join([f"pow(cs[{j}],{s})" for j,s in enumerate(c['nu_f'][:,i]) if s != 0.0])};
  % if c['reversible'][i] == 1.0:
  fpdtype_t q_b = -k_f/K_c * ${"*".join([f"pow(cs[{j}],{s})" for j,s in enumerate(c['nu_b'][:,i]) if s != 0.0])};
  q[${i}] = q_f + q_b;
  % else:
  q[${i}] = q_f;
  % endif
  }
% endfor

  // Output source terms
  src[0] = 0.0;
% for i in range(ndims):
  src[${i+1}] = 0.0;
% endfor
  src[${ndims+1}] = 0.0;
  // Chemical source terms
% for n in range(ns-1):
    // ${c['names'][n]}
<% nu_sum = c['nu_b'][n,:] - c['nu_f'][n,:] %>\
  src[${Yix+n}] = ${MW[n]}*(${"+".join([f"({s}*q[{j}])" for j,s in enumerate(nu_sum) if s != 0.0])});
% endfor

#ifdef DEBUG
  printf("*********************************\n");
  printf("CHEMICAL SOURCE TERMS\n");
% for n in range(ns):
  printf("omega_${c['names'][n]} = %.14f\n", src[${Yix+n}]);
% endfor
  printf("*********************************\n");

#endif

</%pyfr:macro>