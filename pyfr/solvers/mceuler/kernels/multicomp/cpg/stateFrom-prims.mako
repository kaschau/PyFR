<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='stateFrom-prims' params='u, q, qh'>
    ## q is the primitive state, length nvars + 2; indices from mcf.mcix
    ## 0:ns,       ns:ns+ndims, ns+ndims, ns+ndims+1, ns+ndims+2
    ## Y0..Yns-1,  u,v(,w),     rho,      p,          T

    ## qh is the mixture properties, length 4 + ns
    ## 0,     1,  2, 3,     4 + n
    ## gamma, cp, c, rho*e, h_n

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(mcf.ns):
    R += q[${n}] * ${mcf.Ru / mcf[n].MW};
    cp += q[${n}] * ${mcf[n].cp0};
% endfor

    // Compute density
    fpdtype_t rho = q[${pix}] / (R*q[${Tix}]);
    q[${rhoix}] = rho;

    // Species mass
% for n in range(mcf.ns):
    u[${n}] = q[${n}] * rho;
% endfor

    // Compute momentum
% for i in range(ndims):
    u[${i + vix}] = q[${i + vix}] * rho;
% endfor

    // Total energy
    fpdtype_t e = (cp - R) * q[${Tix}];
    u[${Eix}] = rho*(e + 0.5 * ${pyfr.dot('q[{i}]', i=(vix,vix + ndims))});

    // Store gamma, cp, c, rhoe, hs
    qh[0] = cp / (cp - R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0] * R * q[${Tix}]);
    qh[3] = rho * e;

% for n in range(mcf.ns):
    qh[${4 + n}] = q[${Tix}] * ${mcf[n].cp0};
% endfor

% if mcf.debug:
  printf("*********************************\n");
  printf("PRIMS TO CONS\n");
  printf("INPUT STATE\n");
  printf("therm&p = %e\n", q[${pix}]);
% for i in range(ndims):
  printf("therm&v${i + vix} = %e\n", q[${i + vix}]);
% endfor
  printf("therm&T = %e\n", q[${Tix}]);
% for n in range(mcf.ns):
  printf("therm&Y_${mcf.sp_names[n]} = %e\n", q[${n}]);
% endfor

  printf("\nCOMPUTED STATE\n");
  printf("therm&rho = %e\n", q[${rhoix}]);
% for i in range(ndims):
  printf("therm&rhov${i + vix} = %e\n", u[${i + vix}]);
% endfor
  printf("therm&rhoE = %e\n", u[${Eix}]);
% for n in range(mcf.ns):
  printf("therm&rhoY_${mcf.sp_names[n]} = %e\n", u[${n}]);
% endfor
  printf("*********************************\n");
% endif

</%pyfr:macro>