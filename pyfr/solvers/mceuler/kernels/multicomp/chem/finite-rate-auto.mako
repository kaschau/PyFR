<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.chem.net-rate-of-production'/>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='finite_rate_auto' params='t, u, ploc, src'>

  fpdtype_t q[${nvars + 2}];
  fpdtype_t qh[${4 + mcf.ns}];
  ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};

  fpdtype_t rho = q[${rhoix}];
  fpdtype_t rhoinv = 1.0/rho;
  fpdtype_t T = q[${Tix}];

  fpdtype_t tmpSrc[${mcf.ns}];

  % for n in range(mcf.ns):
    src[${n}] = 0.0;
  % endfor

  fpdtype_t tRemaining = 1.0;
  for (int iter = 0; iter < ${max_subs} && tRemaining > ${fpdtype_eps}; iter++){

    ${pyfr.expand('net_rate_of_production', 'q', 'T', 'rho', 'tmpSrc')};

    // Adaptive sub-step size
    fpdtype_t tSubRatio = tRemaining;
    % for n in range(mcf.ns):
<%  nu_sum = mcf.nu_b[n,:] - mcf.nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
    {
      fpdtype_t s = copysign(1.0, tmpSrc[${n}]);
      fpdtype_t headroom = 0.5*((1.0 - s)*q[${n}] + (1.0 + s)*(1.0 - q[${n}]));
      fpdtype_t dtSubMax = rho * fmax(headroom, ${fpdtype_eps}) / fmax(fabs(tmpSrc[${n}]), ${fpdtype_min});
      tSubRatio = fmin(tSubRatio, dtSubMax * ${1.0/dt});
    }
      % endif
    % endfor

    // Compute cp
    fpdtype_t cp = 0.0;
    % for n in range(mcf.ns):
    {
      fpdtype_t cps = ${mcf[n].cp_expr('T')};
      cp += cps*q[${n}];
    }
    % endfor

    fpdtype_t tSub = tSubRatio * ${dt};

    // Temperature and species sub-step
    fpdtype_t dTdt = 0.0;
    % for n in range(mcf.ns):
<%  nu_sum = mcf.nu_b[:,n] - mcf.nu_f[:,n] %>\
      % if max(abs(nu_sum)) > 0.0:
    {
      fpdtype_t hs = ${mcf[n].h_expr('T')};
      dTdt -= hs * tmpSrc[${n}];
      q[${n}] += tmpSrc[${n}] * rhoinv * tSub;
    }
      % endif
      src[${n}] += tmpSrc[${n}] * tSubRatio;
    % endfor

    dTdt /= cp * rho;
    T += dTdt * tSub;

    tRemaining -= tSubRatio;
  }

  % for i in range(ndims):
    src[${i + vix}] = 0.0;
  % endfor
    src[${Eix}] = 0.0;

#ifdef DEBUG
    printf("*********************************\n");
    printf("CHEMICAL SOURCE TERMS (auto)\n");
  % for n in range(mcf.ns):
    printf("chem&omega_${mcf.sp_names[n]} = %e\n", src[${n}]);
  % endfor
    printf("*********************************\n");
#endif

</%pyfr:macro>
