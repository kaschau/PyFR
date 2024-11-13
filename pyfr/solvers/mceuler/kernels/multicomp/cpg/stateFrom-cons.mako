<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='stateFrom-cons' params='u, q, qh'>
    ## q is an array of length nvars + 2
    ## storing all primatives
    ## 0:ns-1,    ns:ns+ndims, ns+ndims+1, ns+ndims+2, nvars + 2
    ## Y0...Ynsp, u,v(,w),  rho          , p,          T

    ## qh stores mixture thermodynamic properties
    ## 0,  1,     2, 3, 4..4 + ns
    ## cp, gamma, c, e, hi1..hins

    // Compute rho
    fpdtype_t rho = ${" + ".join([f"u[{n}]" for n in range(ns)])};
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${Eix}];

    // Compute species mass fraction
% for n in range(ns):
    q[${n}] = u[${n}]*invrho;
% endfor

    // Compute velocities
% for i in range(ndims):
    q[${i + vix}] = u[${i + vix}]*invrho;
% endfor

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    R += q[${n}]*${1.0/c['MW'][n]};
    cp += q[${n}]*${c['cp0'][n]};
% endfor
    R *= ${c['Ru']};

    // Internal energy (per mass)
    fpdtype_t e = (rhoE - 0.5*rho*${pyfr.dot('q[{i}]', i=(vix,vix + ndims))})*invrho;

    // Equilibrium T, p
    q[${rhoix}] = rho;
    q[${Tix}] = e / (cp - R);
    q[${pix}] = rho * R * q[${Tix}];

    // Mixture gamma, cp
    qh[0] = cp / (cp - R);
    qh[1] = cp;
    // Mixture speed of sound
    qh[2] = sqrt(qh[0] * R * q[${Tix}]);
    // internal energy
    qh[3] = rho*e;

    // Store species enthalpy (per mass)
% for n in range(ns):
    qh[${4 + n}] = q[${Tix}]*${c['cp0'][n]};
% endfor

#ifdef DEBUG
  printf("*********************************\n");
  printf("CPG THERMODYNAMIC PROPERTIES\n");
  printf("INPUT STATE\n");
% for n in range(ns):
  printf("therm&rhoY_${c['names'][n]} = %f\n", u[${n}]);
% endfor
% for i in range(ndims):
    printf("therm&rhou = %f\n", u[${vix + i}]);
% endfor
  printf("therm&rhoE = %f\n", u[${Eix}]);

  printf("\nCOMPUTED STATE\n");
  printf("therm&rho = %f\n", q[${rhoix}]);
  printf("therm&p = %f\n", q[${pix}]);
  printf("therm&T = %f\n", q[${Tix}]);
% for n in range(ns):
  printf("therm&Y_${c['names'][n]} = %f\n", q[${n}]);
% endfor
% for i in range(ndims):
    printf("therm&u = %f\n", q[${vix + i}]);
% endfor

  printf("\nCOMPUTED PROPERTIES\n");
  printf("therm&R = %f\n", R);
  printf("therm&gamma = %f\n", qh[0]);
  printf("therm&cp = %f\n", qh[1]);
  printf("therm&c = %f\n", qh[2]);
  printf("therm&rhoe = %f\n", qh[3]);
% for n in range(ns):
  printf("therm&h_${c['names'][n]} = %f\n", qh[${4 + n}]);
% endfor
  printf("*********************************\n");
#endif
</%pyfr:macro>