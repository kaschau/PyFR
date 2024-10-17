<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

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
<% nr = c['Ea_f'].shape[0] %>\
<% MW = c['MW'] %>\
<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% div = [1.0, 2.0, 3.0, 4.0, 5.0] %>\

  // Compute thermodynamic properties
  fpdtype_t q[${nvars + 2}];
  fpdtype_t qh[${4 + ns}];
  ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};

  fpdtype_t rho = q[${rhoix}];
  fpdtype_t tSub = ${dt} / ${nsub_steps};

  fpdtype_t T = q[${Tix}];

  for(int nSub = 0; nSub < ${nsub_steps}; nSub++){

  // Concentrations
  fpdtype_t cs[${ns}];
% for n in range(ns):
  cs[${n}] = rho*q[${n}]*${1.0/c['MW'][n]};
% endfor

  // Gibbs energy
  fpdtype_t gbs[${ns}];
  fpdtype_t logT = log(T);
  fpdtype_t Tinv = 1.0/T;
  fpdtype_t prefRuT = ${101325.0/c['Ru']}*Tinv;
  fpdtype_t cp = 0.0;
  {
% for n in range(ns):
      // ${c['names'][n]} Properties
      if (T < ${N7[n,0]})
      {
<% m = 8 %>\
            fpdtype_t cps = ${'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n])+')'*4};
            fpdtype_t hi = ${f'+ T*('.join(str(c) for c in N7[n,m:m+5]/div)+')'*4} + ${N7[n, m + 5]}*Tinv;
            fpdtype_t scs = ${N7[n, m + 0]} * logT +
                            T*(${f'+ T*('.join(str(c) for c in N7[n,m+1:m+5]/div[0:-1])+')'*3}) + ${N7[n, m + 6]};
            gbs[${n}] = hi - scs;
            qh[${4 + n}] = hi;
            cp += cps*q[${n}];
        }else
        {
<% m = 1 %>\
            fpdtype_t cps = ${'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n])+')'*4};
            fpdtype_t hi = ${f'+ T*('.join(str(c) for c in N7[n,m:m+5]/div)+')'*4} + ${N7[n, m + 5]}*Tinv;
            fpdtype_t scs = ${N7[n, m + 0]} * logT +
                            T*(${f'+ T*('.join(str(c) for c in N7[n,m+1:m+5]/div[0:-1])+')'*3}) + ${N7[n, m + 6]};
            gbs[${n}] = hi - scs;
            qh[${4 + n}] = hi;
            cp += cps*q[${n}];
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
  fpdtype_t rp[${nr}];

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
  % if c['r_type'][i] == "Arrhenius Custom Order":
  fpdtype_t rp_f = k_f * ${"*".join([f"pow(cs[{j}],{s})" for j,s in enumerate(c['orders'][i]) if s != 0.0])};
  % else:
  fpdtype_t rp_f = k_f * ${"*".join([f"pow(cs[{j}],{s})" for j,s in enumerate(c['nu_f'][:,i]) if s != 0.0])};
  % endif
  % if c['reversible'][i] == 1.0:
  fpdtype_t rp_b = -k_f/K_c * ${"*".join([f"pow(cs[{j}],{s})" for j,s in enumerate(c['nu_b'][:,i]) if s != 0.0])};
  rp[${i}] = rp_f + rp_b;
  % else:
  rp[${i}] = rp_f;
  % endif
  }
% endfor

  fpdtype_t dTdt = 0.0;
  fpdtype_t tempsum = 0.0;
% for n in range(ns):
  {
<% nu_sum = c['nu_b'][n,:] - c['nu_f'][n,:] %>\
% if max(abs(nu_sum)) > 0:
    fpdtype_t dYdt = ${MW[n]}*(${"+".join([f"({s}*rp[{j}])" for j,s in enumerate(nu_sum) if s != 0.0])});
% else:
    fpdtype_t dYdt = 0.0;
% endif
    dTdt -= qh[${4 + n}] * dYdt;
    q[${n}] += dYdt / rho * tSub;
    q[${n}] = fmax(0.0, q[${n}]);
    tempsum += q[${n}];
  }
% endfor
// Normalize
% for n in range(ns):
    q[${n}] /= tempsum;
% endfor

  dTdt /= cp*rho;
  T += dTdt*tSub;
  }

  // Chemical source terms
  // Reconstruct d(rhoY)/dt based on where we ended up
% for n in range(ns):
  // ${c['names'][n]}
  src[${n}] = (q[${n}] * rho - u[${n}]) / ${dt};
% endfor

% for i in range(ndims):
  src[${i + vix}] = 0.0;
% endfor

  src[${Eix}] = 0.0;


#ifdef DEBUG
  printf("*********************************\n");
  printf("CHEMICAL SOURCE TERMS\n");
% for n in range(ns):
  printf("chem&omega_${c['names'][n]} = %.14f\n", src[${n}]);
% endfor
  printf("*********************************\n");
#endif

</%pyfr:macro>