<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.tpg.T_iter' />

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='stateFrom-cons' params='u, q, qh'>

    // Compute rho
    fpdtype_t rho = ${" + ".join([f"u[{n}]" for n in range(mcf.ns)])};
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${Eix}];

    // Mass fractions
% for n in range(mcf.ns):
    q[${n}] = u[${n}]*invrho;
% endfor

    // Compute velocities
% for i in range(ndims):
    q[${i + vix}] = u[${i + vix}]*invrho;
% endfor

    // Compute mixture gas constant
    fpdtype_t R = 0.0;
% for n in range(mcf.ns):
    R += q[${n}]*${mcf.Ru / mcf[n].MW};
% endfor

    // Internal energy (per mass)
    fpdtype_t e = rhoE*invrho - 0.5*${pyfr.dot('q[{i}]', i=(vix,vix + ndims))};

    // Iterate on T
    fpdtype_t cp;
    fpdtype_t T;
    ${pyfr.expand(f'T_iter_{mcf.T_iter_method}', 'e', 'cp', 'R', 'T', 'q', 'qh')};

    // Equilibrium T, p
    q[${rhoix}] = rho;
    q[${pix}] = rho * R * T;
    q[${Tix}] = T;

    // Mixture gamma, cp
    qh[0] = cp / (cp - R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0] * R * T);
    qh[3] = rho*e;

#ifdef DEBUG
  printf("*********************************\n");
  printf("TPG THERMODYNAMIC PROPERTIES\n");
  printf("INPUT STATE\n");
  printf("therm&rho = %e\n", q[${rhoix}]);
  printf("therm&e = %e\n", qh[3]/rho);
% for n in range(mcf.ns):
  printf("therm&rhoY_${mcf.sp_names[n]} = %e\n", u[${n}]);
% endfor

  printf("\nCOMPUTED STATE\n");
  printf("therm&p = %e\n", q[${pix}]);
  printf("therm&T = %e\n", q[${Tix}]);
% for n in range(mcf.ns):
  printf("therm&Y_${mcf.sp_names[n]} = %e\n", q[${n}]);
% endfor

  printf("\nCOMPUTED PROPERTIES\n");
  printf("therm&R = %e\n", R);
  printf("therm&gamma = %e\n", qh[0]);
  printf("therm&cp = %e\n", qh[1]);
  printf("therm&c = %e\n", qh[2]);
% for n in range(mcf.ns):
  printf("therm&h_${mcf.sp_names[n]} = %e\n", qh[${4 + n}]);
% endfor
  printf("*********************************\n");
#endif
</%pyfr:macro>
