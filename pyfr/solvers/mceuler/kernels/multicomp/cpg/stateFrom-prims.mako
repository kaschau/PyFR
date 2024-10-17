<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='stateFrom-prims' params='u, q, qh'>
    ## q is an array of length nvars + 2
    ## storing all primatives
    ## 0:ns-1,    ns:ns+ndims, ns+ndims+1, ns+ndims+2, nvars + 2
    ## Y0...Ynsp, u,v(,w),  rho          , p,          T

    ## qh stores mixture thermodynamic properties
    ## 0,  1,     2, 3, 4..4 + ns
    ## cp, gamma, c, e, hi1..hins

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    R += q[${n}]*${1.0/c['MW'][n]};
    cp += q[${n}]*${c['cp0'][n]};
% endfor
    R *= ${c['Ru']};

    // Compute density
    fpdtype_t rho = q[${pix}]/(R*q[${Tix}]);
    q[${rhoix}] = rho;

    // Species mass
% for n in range(ns):
    u[${n}] = q[${n}]*rho;
% endfor

    // Compute momentum
% for i in range(ndims):
    u[${i + vix}] = q[${i + vix}]*rho;
% endfor

    // Total energy
    fpdtype_t e = (cp-R) * q[${Tix}];
    u[${Eix}] = rho*(e + 0.5*${pyfr.dot('q[{i}]', i=(vix,vix + ndims))});

    // Store gamma, cp, c, rhoe, hs
    qh[0] = cp/(cp-R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0]*R*q[${Tix}]);
    qh[3] = rhoe
% for n in range(ns):
    qh[${4 + n}] = q[${Tix}]*${c['cp0'][n]};
% endfor

#ifdef DEBUG
  printf("*********************************\n");
  printf("PRIMS TO CONS");
  printf("INPUT STATE\n");
  printf("therm&p = %.14f\n", q[${pix}]);
% for i in range(ndims):
  printf("therm&v${i + vix} = %.14f\n", q[${i + vix}]);
% endfor
  printf("therm&T = %.14f\n", q[${Tix}]);
% for n in range(ns):
  printf("therm&Y_${c['names'][n]} = %.14f\n", q[${n}]);
% endfor

  printf("\nCOMPUTED STATE\n");
  printf("therm&rho = %.14f\n", q[${rhoix}]);
% for i in range(ndims):
  printf("therm&rhov${i + vix} = %.14f\n", u[${i + vix}]);
% endfor
  printf("therm&rhoE = %.14f\n", u[${Eix}]);
% for n in range(ns):
  printf("therm&rhoY_${c['names'][n]} = %.14f\n", u[${n}]);
% endfor
  printf("*********************************\n");
#endif

</%pyfr:macro>