<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.tpg.T_iter'/>

<%pyfr:macro name='mixture_state' params='u, q, qh'>
    ## q is an array of length nvars + 1
    ## storing all primatives
    ## 0, 1:ndims, ndims+1, ndims+2 ... nvars
    ## p, u,v(,w),   T    ,   Y0    ...  Yns

    ## qh stores mixture thermodynamic properties
    ## 0,  1,     2, 3..3+ns
    ## cp, gamma, c, hi1..hins
    fpdtype_t rho = u[0];
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${ndims + 1}];

    // Compute velocities
% for i in range(ndims):
    q[${i+1}] = u[${i+1}]*invrho;
% endfor

## First species index
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>

    q[${nvars}] = 1.0;
% for n in range(ns-1):
    q[${Yix+n}] = u[${Yix+n}]*invrho;
    q[${nvars}] -= q[${Yix+n}];
% endfor

    // Compute mixture properties
    fpdtype_t R = 0.0;
% for n in range(ns):
    R += q[${Yix+n}]*${1.0/c['MW'][n]};
% endfor
    R *= ${c['Ru']};

    // Internal energy (per mass)
    fpdtype_t e = rhoE*invrho - 0.5*${pyfr.dot('q[{i}]', i=(1,ndims+1))};

    // Iterate on T
    fpdtype_t h;
    fpdtype_t cp;
    fpdtype_t T = 300.0; // Initial guess
    ${pyfr.expand('T_iter', 'e', 'h', 'cp', 'R', 'T', 'q', 'qh')};

    // Equilibrium T, p
    q[0] = rho * R * T;
    q[${ndims + 1}] = T;

    // Mixture gamma, cp
    qh[0] = cp / (cp - R);
    qh[1] = cp;

    // Mixture speed of sound
    qh[2] = sqrt(qh[0] * R * q[${ndims + 1}]);

    // Store species enthalpy (per mass)
    // ^ done in T_iter

#ifdef DEBUG
  printf("*********************************\n");
  printf("THERMODYNAMIC PROPERTIES\n");
  printf("INPUT STATE\n");
  printf("rho = %.14f\n", u[0]);
  printf("e = %.14f\n", e);
% for n in range(ns-1):
  printf("rhoY_${c['names'][n]} = %.14f\n", u[${Yix+n}]);
% endfor

  printf("\nCOMPUTED STATE\n");
  printf("p = %.14f\n", q[0]);
  printf("T = %.14f\n", q[${ndims+1}]);
% for n in range(ns):
  printf("Y_${c['names'][n]} = %.14f\n", q[${Yix+n}]);
% endfor

  printf("\nCOMPUTED PROPERTIES\n");
  printf("gamma = %.14f\n", qh[0]);
  printf("cp = %.14f\n", qh[1]);
  printf("c = %.14f\n", qh[2]);
% for n in range(ns):
  printf("h_${c['names'][n]} = %.14f\n", qh[${3+n}]);
% endfor
  printf("*********************************\n");

#endif
</%pyfr:macro>