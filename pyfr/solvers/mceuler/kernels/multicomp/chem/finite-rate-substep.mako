<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.chem.net-rate-of-production'/>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>
<% tSub = dt / float(sub_steps) %>

<%pyfr:macro name='finite_rate_substep' params='t, u, ploc, src'>

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

  for(int nSub = 0; nSub < ${sub_steps}; nSub++){
    fpdtype_t ud[${mcf.ns}];
    % for n in range(mcf.ns):
    ud[${n}] = rho*q[${n}];
    % endfor
    ${pyfr.expand('net_rate_of_production', 'ud', 'T', 'tmpSrc')};

    // Compute cp
    fpdtype_t cp = 0.0;
    % for n in range(mcf.ns):
    {
      fpdtype_t cps = ${mcf[n].cp_expr('T')};
      cp += cps*q[${n}];
    }
    % endfor

    // Temperature sub-step
    fpdtype_t dTdt = 0.0;
    % for n in range(mcf.ns):
<%  nu_sum = mcf.nu_b[n,:] - mcf.nu_f[n,:] %>\
      % if max(abs(nu_sum)) > 0.0:
    {
      fpdtype_t hs = ${mcf[n].h_expr('T')};
      dTdt -= hs * tmpSrc[${n}];
    }
      % endif
      src[${n}] += tmpSrc[${n}]*${tSub/dt};
    % endfor
    dTdt /= cp * rho;
    T += dTdt * ${tSub};

    // Species sub-step (clamp + normalize)
    fpdtype_t Y_sum = 0.0;
    % for n in range(mcf.ns):
<%  nu_sum = mcf.nu_b[:,n] - mcf.nu_f[:,n] %>\
      % if max(abs(nu_sum)) > 0.0:
      q[${n}] = fmin(1.0, fmax(0.0, q[${n}] + tmpSrc[${n}] * rhoinv * ${tSub}));
      % endif
      Y_sum += q[${n}];
    % endfor
    fpdtype_t Y_suminv = 1.0/Y_sum;
    % for n in range(mcf.ns):
      q[${n}] *= Y_suminv;
    % endfor
  }

  % for i in range(ndims):
    src[${i + vix}] = 0.0;
  % endfor
    src[${Eix}] = 0.0;

#ifdef DEBUG
    printf("*********************************\n");
    printf("CHEMICAL SOURCE TERMS (substep)\n");
  % for n in range(mcf.ns):
    printf("chem&omega_${mcf.sp_names[n]} = %e\n", src[${n}]);
  % endfor
    printf("*********************************\n");
#endif

</%pyfr:macro>
