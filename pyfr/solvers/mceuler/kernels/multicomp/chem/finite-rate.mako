<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.chem.net-rate-of-production'/>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='finite_rate' params='t, u, ploc, src'>

  fpdtype_t q[${nvars + 2}];
  fpdtype_t qh[${4 + mcf.ns}];
  ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};

  fpdtype_t rho = q[${rhoix}];
  fpdtype_t T = q[${Tix}];

  ${pyfr.expand('net_rate_of_production', 'q', 'T', 'rho', 'src')};

  % for i in range(ndims):
    src[${i + vix}] = 0.0;
  % endfor
    src[${Eix}] = 0.0;

#ifdef DEBUG
    printf("*********************************\n");
    printf("CHEMICAL SOURCE TERMS\n");
  % for n in range(mcf.ns):
    printf("chem&omega_${mcf.sp_names[n]} = %e\n", src[${n}]);
  % endfor
    printf("*********************************\n");
#endif

</%pyfr:macro>
