<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='stateFrom-prims' params='u, q, qh'>
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>

    ## q is an array of length nvars + 1
    ## storing all primatives
    ## 0, 1:ndims, ndims+1, ndims+2 ... nvars
    ## p, u,v(,w),   T    ,   Y0    ...  Yns

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    R += q[${Yix+n}]*${1.0/c['MW'][n]};
    cp += q[${Yix+n}]*${c['cp0'][n]};
% endfor
    R *= ${c['Ru']};

    // Compute density
    fpdtype_t rho = q[0]/(R*q[${ndims+1}]);
    u[0] = rho;

    // Compute momentum
% for i in range(ndims):
    u[${i+1}] = q[${i+1}]*rho;
% endfor

    // Total energy
    fpdtype_t e = (cp-R) * q[${ndims+1}];
    u[${ndims+1}] = rho*(e + 0.5*${pyfr.dot('q[{i}]', i=(1,ndims+1))});

    // Species mass
% for n in range(ns-1):
    u[${Yix+n}] = q[${Yix+n}]*rho;
% endfor

    // Store gamma, cp, c, rhoe, hs
    qh[0] = cp/(cp-R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0]*R*q[${ndims+1}]);
    qh[3] = rhoe
% for n in range(ns):
    qh[${4+n}] = q[${ndims+1}]*${c['cp0'][n]};
% endfor

#ifdef DEBUG
  printf("*********************************\n");
  printf("PRIMS TO CONS");
  printf("INPUT STATE\n");
  printf("therm&p = %.14f\n", q[0]);
% for i in range(ndims):
  printf("therm&v${i+1} = %.14f\n", q[${i+1}]);
% endfor
  printf("therm&T = %.14f\n", q[${ndims+1}]);
% for n in range(ns-1):
  printf("therm&Y_${c['names'][n]} = %.14f\n", q[${Yix+n}]);
% endfor

  printf("\nCOMPUTED STATE\n");
  printf("therm&rho = %.14f\n", u[0]);
% for i in range(ndims):
  printf("therm&rhov${i+1} = %.14f\n", u[${i+1}]);
% endfor
  printf("therm&rhoE = %.14f\n", u[${ndims+1}]);
% for n in range(ns-1):
  printf("therm&rhoY_${c['names'][n]} = %.14f\n", u[${Yix+n}]);
% endfor
  printf("*********************************\n");
#endif

</%pyfr:macro>