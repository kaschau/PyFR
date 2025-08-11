<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.chem.net-rate-of-production'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\
<% MW = c['MW'] %>\
<% Ru = c['Ru'] %>\
<% nu_f = c['nu_f'] %>\
<% nu_b = c['nu_b'] %>\
<% fast_props = 'fast_coeff' in c %>\
% if fast_props:
<% fast_coeff = c['fast_coeff'] %>\
% else:
<% T_cutoff = c['T_cutoff'] %>\
<% NASA7_Thigh = c['NASA7_Thigh'] %>\
<% NASA7_Tlow = c['NASA7_Tlow'] %>\
% endif
<% tSub = dt / float(sub_steps) %>\

<%pyfr:macro name='finite_rate_substep' params='t, u, ploc, src'>

  // Compute thermodynamic properties
  fpdtype_t q[${nvars + 2}];
  fpdtype_t qh[${4 + ns}];
  ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};

  fpdtype_t rho = q[${rhoix}];
  fpdtype_t rhoinv = 1.0/rho;
  fpdtype_t T = q[${Tix}];

  fpdtype_t tmpSrc[${ns}];

  // Start source at zero
  % for n in range(ns):
    src[${n}] = 0.0;
  % endfor

  for(int nSub = 0; nSub < ${sub_steps}; nSub++){
    ${pyfr.expand('net_rate_of_production', 'q', 'T', 'rho', 'tmpSrc')};

    // Compute cp
    fpdtype_t cp = 0.0;
    % for n in range(ns):
    {
      % if fast_props:
      {
        fpdtype_t cps = ${pyfr.nasa_cps(fast_coeff[n], Ru, MW[n])};
        cp += cps*q[${n}];
      }
      % else:
      if (T < ${T_cutoff[n]})
      {
        fpdtype_t cps = ${pyfr.nasa_cps(NASA7_Tlow[n], Ru, MW[n])};
        cp += cps*q[${n}];
      }else{
        fpdtype_t cps = ${pyfr.nasa_cps(NASA7_Thigh[n], Ru, MW[n])};
        cp += cps*q[${n}];
      }
      % endif
    }
    % endfor

    // Take sub step in time for temperature
    fpdtype_t dTdt = 0.0;
    fpdtype_t Tinv = 1.0/T;
    % for n in range(ns):
    {
      <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
        % if fast_props:
          fpdtype_t hi = ${pyfr.nasa_hi(fast_coeff[n])};
        % else:
          fpdtype_t hi;
          if (T < ${T_cutoff[n]})
          {
            hi = ${pyfr.nasa_hi(NASA7_Tlow[n])};
          }else
          {
            hi = ${pyfr.nasa_hi(NASA7_Thigh[n])};
          }
        % endif
        dTdt -= hi * tmpSrc[${n}];
      % endif

      // Accumulate source term
      src[${n}] += tmpSrc[${n}]*${tSub/dt};
    }
    % endfor
    dTdt /= cp * rho;
    T += dTdt * ${tSub};

    // Take sub step in time for species
    fpdtype_t Y_sum = 0.0;
    % for n in range(ns):
    {
      <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
        q[${n}] = fmin(1.0, fmax(0.0, q[${n}] + tmpSrc[${n}] * rhoinv * ${tSub}));
      % endif
        Y_sum += q[${n}];
    }
    % endfor
    // Normalize
    fpdtype_t Y_suminv = 1.0/Y_sum;
    % for n in range(ns):
      q[${n}] *= Y_suminv;
    % endfor
  }

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