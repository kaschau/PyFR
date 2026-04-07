<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%namespace module='pyfr.multicomp.makoutil' name='mc'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.chem.net-rate-of-production'/>

<% ns, vix, Eix, rhoix, pix, Tix = mc.thermix(c['ns'], ndims) %>\
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
% endif\

<%pyfr:macro name='finite_rate_auto' params='t, u, ploc, src'>

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

  fpdtype_t tProgress = 0.0;  // Dimensionless progress [0-1]
  for (int iter = 0; iter < ${max_subs} && tProgress < 1.0; iter++){

    ${pyfr.expand('net_rate_of_production', 'q', 'T', 'rho', 'tmpSrc')};

    // Find largest possible sub step ratio (dimensionless)
    fpdtype_t tSubRatio = 1.0 - tProgress;
    % for n in range(ns):
    {
      <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
      ## g.t.zero and l.t. one
      // HACK underflows in tmpSrc fail the sign check with ffast-math
      // resulting in tSub=inf and such. Need better solution.
      if (abs(tmpSrc[${n}]) > ${fpdtype_eps}){
        fpdtype_t dtSubMax = (tmpSrc[${n}] < 0.0) ? -rho*q[${n}]/tmpSrc[${n}]
                                                  : rho*(1.0-q[${n}])/tmpSrc[${n}];
        tSubRatio = fmin(tSubRatio, dtSubMax * ${1.0/dt});
      }
      % endif
    }
    % endfor

    // Compute cp
    fpdtype_t cp = 0.0;
    % for n in range(ns):
    {
      % if fast_props:
      {
        fpdtype_t cps = ${mc.nasa_cps(fast_coeff[n], Ru, MW[n])};
        cp += cps*q[${n}];
      }
      % else:
      if (T < ${T_cutoff[n]})
      {
        fpdtype_t cps = ${mc.nasa_cps(NASA7_Tlow[n], Ru, MW[n])};
        cp += cps*q[${n}];
      }else{
        fpdtype_t cps = ${mc.nasa_cps(NASA7_Thigh[n], Ru, MW[n])};
        cp += cps*q[${n}];
      }
      % endif
    }
    % endfor

    // Convert ratio to actual time step
    fpdtype_t tSub = tSubRatio * ${dt};

    // Take the sub-step
    fpdtype_t dTdt = 0.0;
    fpdtype_t Tinv = 1.0/T;
    % for n in range(ns):
    {
      <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
        % if fast_props:
          fpdtype_t hi = ${mc.nasa_hi(fast_coeff[n])};
        % else:
          fpdtype_t hi;
          if (T < ${T_cutoff[n]})
          {
            hi = ${mc.nasa_hi(NASA7_Tlow[n])};
          }else
          {
            hi = ${mc.nasa_hi(NASA7_Thigh[n])};
          }
        % endif
        dTdt -= hi * tmpSrc[${n}];
      % endif

      // Take sub step in time for species
      <% nu_sum = nu_b[n,:] - nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
        q[${n}] += tmpSrc[${n}] * rhoinv * tSub;
      % endif

      // Accumulate source term
      src[${n}] += tmpSrc[${n}] * tSubRatio;
    }
    % endfor

    // Take sub step in time for temperature
    dTdt /= cp * rho;
    T += dTdt * tSub;

    tProgress += tSubRatio;
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