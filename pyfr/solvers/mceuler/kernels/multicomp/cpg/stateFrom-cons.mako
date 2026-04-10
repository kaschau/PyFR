<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='stateFrom-cons' params='u, q, qh'>
    ## q is an array of length nvars + 2
    ## storing all primitive s
    ## 0:ns-1,    ns:ns+ndims, ns+ndims+1, ns+ndims+2, nvars + 2
    ## Y0...Ynsp, u,v(,w),  rho          , p,          T

    ## qh stores mixture thermodynamic properties
    ## 0,  1,     2, 3, 4..4 + ns
    ## cp, gamma, c, e, hi1..hins

    // Compute rho
    fpdtype_t rho = ${" + ".join([f"u[{n}]" for n in range(mcf.ns)])};
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${Eix}];

    // Compute species mass fraction
% for n in range(mcf.ns):
    q[${n}] = u[${n}]*invrho;
% endfor

    // Compute velocities
% for i in range(ndims):
    q[${i + vix}] = u[${i + vix}]*invrho;
% endfor

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(mcf.ns):
    R += q[${n}]*${mcf.Ru / mcf[n].MW};
    cp += q[${n}]*${mcf[n].cp0};
% endfor

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
% for n in range(mcf.ns):
    qh[${4 + n}] = q[${Tix}]*${mcf[n].cp0};
% endfor

#ifdef DEBUG
  printf("*********************************\n");
  printf("CPG THERMODYNAMIC PROPERTIES\n");
  printf("INPUT STATE\n");
% for n in range(mcf.ns):
  printf("therm&rhoY_${mcf.sp_names[n]} = %e\n", u[${n}]);
% endfor
% for i in range(ndims):
    printf("therm&rhou = %e\n", u[${vix + i}]);
% endfor
  printf("therm&rhoE = %e\n", u[${Eix}]);

  printf("\nCOMPUTED STATE\n");
  printf("therm&rho = %e\n", q[${rhoix}]);
  printf("therm&p = %e\n", q[${pix}]);
  printf("therm&T = %e\n", q[${Tix}]);
% for n in range(mcf.ns):
  printf("therm&Y_${mcf.sp_names[n]} = %e\n", q[${n}]);
% endfor
% for i in range(ndims):
    printf("therm&u = %e\n", q[${vix + i}]);
% endfor

  printf("\nCOMPUTED PROPERTIES\n");
  printf("therm&R = %e\n", R);
  printf("therm&gamma = %e\n", qh[0]);
  printf("therm&cp = %e\n", qh[1]);
  printf("therm&c = %e\n", qh[2]);
  printf("therm&rhoe = %e\n", qh[3]);
% for n in range(mcf.ns):
  printf("therm&h_${mcf.sp_names[n]} = %e\n", qh[${4 + n}]);
% endfor
  printf("*********************************\n");
#endif
</%pyfr:macro>