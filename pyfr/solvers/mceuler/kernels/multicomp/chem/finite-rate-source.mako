<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<% import math %>\

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\

<% nr = c['Ea_f'].shape[0] %>\
<% MW = c['MW'] %>\
<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% fast_props = N7.shape[1] == 7 %>\
<% reconstruct = nsub_steps > 1 %>\
<% tSub = dt / float(nsub_steps) %>\

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
    prefRuTinv*exp(dG)
  % elif nusum == -1.0:
    exp(dG)*prefRuT
  % elif nusum.is_integer():
    % if nusum > 0.0:
      exp(dG)*${pyfr.intpow("prefRuT", nusum)}
    % else:
      exp(dG)*${pyfr.intpow("prefRuT", -(nusum))}
    % endif:
  % else:
    pow(prefRuT,-(${nusum}))*exp(dG)
  %   endif
% else:
    exp(dG)
% endif
</%def>\

<%pyfr:macro name='finite_rate_source' params='t, u, ploc, src'>

  // Compute thermodynamic properties
  fpdtype_t q[${nvars + 2}];
  fpdtype_t qh[${4 + ns}];
  ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};

  fpdtype_t rho = q[${rhoix}];

  fpdtype_t T = q[${Tix}];

  for(int nSub = 0; nSub < ${nsub_steps}; nSub++){

  // Concentrations
  fpdtype_t cs[${ns}];
  % for n in range(ns):
    cs[${n}] = rho*q[${n}]*${1.0/c['MW'][n]};
  % endfor

  // Gibbs energy
  fpdtype_t gbs[${ns}];
  double logT = log(T);
  double Tinv = 1.0/T;
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

  fpdtype_t rp[${nr}];
  % for i in range(nr):
  <% alpha = c['fall_coeffs'][i][0]%>\
  <% Tsss = c['fall_coeffs'][i][1]%>\
  <% Ts = c['fall_coeffs'][i][2]%>\
  <% Tss = c['fall_coeffs'][i][3]%>\
  <% nu_sum = nu_b[:,i] - nu_f[:,i] %>\

  // Reaction ${i} - ${c['r_type'][i]}
  {
  double k_f = ${rateConst(A_f[i], m_f[i], Ea_f[i])};
  % if sum(c['aij'][i]) > 0.0:
  // Three body reaction
  fpdtype_t cTBC = ${"+".join([f"({eff}*cs[{j}])" for j,eff in enumerate(c['aij'][i]) if eff != 0.0])};
  % endif
  % if c['r_type'][i] == 'three-body-Arrhenius':
    k_f *= cTBC;
  % elif c['r_type'][i] == 'falloff-Lindemann':
    // Lindemann Reaction
    fpdtype_t Pr = cTBC*${rateConst(A_o[i]/A_f[i], m_o[i]-m_f[i], Ea_o[i]-Ea_f[i])}; // <- ratio k0/k_f
    fpdtype_t pmod = Pr/(1.0 + Pr);
    k_f *= pmod;
  % elif c['r_type'][i] == 'falloff-Troe':
    // Troe Reaction
    % if Tss == 0: #Three Parameter Troe form
      fpdtype_t log10Fcent = log10((${1.0 - alpha})*exp(-T*${1.0/Tsss}) + ${alpha}*exp(-T*${1.0/Ts}));
    % else: # Four Parameter Troe form
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
  % if c['r_type'][i] == "Arrhenius Custom Order":
    rp[${i}] = k_f * ${"*".join([pyfr.intpow(f"cs[{j}]",s) for j,s in enumerate(c['orders'][i]) if float(s) != 0.0])};
  % else:
    rp[${i}] = k_f * ${"*".join([pyfr.intpow(f"cs[{j}]",s) for j,s in enumerate(nu_f[:,i]) if float(s) != 0.0])};
  % endif

  % if c['reversible'][i] == 1.0:
    double dG = ${"+".join([f"({s}*gbs[{i}])" for i,s in enumerate(nu_sum) if s != 0.0])};
    double K_cinv = ${Kcinv(sum(nu_sum))};
    rp[${i}] -= k_f*K_cinv * ${"*".join([pyfr.intpow(f"cs[{j}]",s) for j,s in enumerate(nu_b[:,i]) if float(s) != 0.0])};
  % endif
  }
% endfor ##// End reaction loop

% if reconstruct:
  fpdtype_t cp = 0.0;
  % for n in range(ns):
  {
    % if fast_props:
    {
      fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 0)};
      cp += cps*q[${n}];
    }
    % else:
    if (T < ${N7[n,0]}){
      fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 8)};
      cp += cps*q[${n}];
    }else{
      fpdtype_t cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 1)};
      cp += cps*q[${n}];
    }
    % endif
  }
  % endfor
  // Take sub step in time
  fpdtype_t dTdt = 0.0;
  fpdtype_t tempsum = 0.0;
  fpdtype_t rhoinv = 1.0/rho;
  % for n in range(ns):
  {
    <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
    % if max(abs(nu_sum)) > 0.0:
      fpdtype_t dYdt = ${MW[n]}*(${"+".join([f"({s}*rp[{j}])" for j,s in enumerate(nu_sum) if s != 0.0])});
      % if fast_props:
        fpdtype_t hi = ${pyfr.nasa_hi(N7[n,:], 0)};
      % else:
        fpdtype_t hi;
        if (T < ${N7[n,0]})
        {
          hi = ${pyfr.nasa_hi(N7[n,:], 8)};
        }else
        {
          hi = ${pyfr.nasa_hi(N7[n,:], 1)};
        }
      % endif
      dTdt -= hi * dYdt;
      q[${n}] += dYdt *rhoinv * ${tSub};
      q[${n}] = fmax(0.0, q[${n}]);
    % endif
    tempsum += q[${n}];
  }
  % endfor
  // Normalize
  % for n in range(ns):
    q[${n}] /= tempsum;
  % endfor
  dTdt /= cp * rho;
  T += dTdt * ${tSub};

% else: ## Not reconstructing

  // Chemical source terms
  // Just set the source term
  % for n in range(ns):
    // ${c['names'][n]}
    <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
    % if max(abs(nu_sum)) > 0.0:
      src[${n}] = ${MW[n]}*(${"+".join([f"({s}*rp[{j}])" for j,s in enumerate(nu_sum) if s != 0.0])});
    % else:
      src[${n}] = 0.0;
    % endif
  % endfor

% endif
}

% if reconstruct:
  // Chemical source terms
  // Reconstruct d(rhoY)/dt based on where we ended up
  % for n in range(ns):
    // ${c['names'][n]}
    <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
    % if max(abs(nu_sum)) > 0.0:
      % if reconstruct:
        src[${n}] = (q[${n}] * rho - u[${n}]) / ${dt};
      % endif
    % else:
      src[${n}] = 0.0;
    % endif
  % endfor
% endif

// Set non chemical terms to zero
% for i in range(ndims):
  src[${i + vix}] = 0.0;
% endfor

  src[${Eix}] = 0.0;


#ifdef DEBUG
  printf("*********************************\n");
  printf("CHEMICAL SOURCE TERMS\n");
% for n in range(ns):
  printf("chem&omega_${c['names'][n]} = %e\n", src[${n}]);
% endfor
  printf("*********************************\n");
#endif

</%pyfr:macro>