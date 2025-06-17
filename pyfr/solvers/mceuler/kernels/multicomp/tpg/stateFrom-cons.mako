<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.tpg.T_iter' />

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = N7.shape[1] == 7 %>\

<%pyfr:macro name='stateFrom-cons' params='u, q, qh'>

    ## q is an array of length nvars + 2
    ## storing all primitive s
    ## 0:ns-1,    ns:ns+ndims, ns+ndims+1, ns+ndims+2, nvars + 2
    ## Y0...Ynsp, u,v(,w),  rho          , p,          T

    ## qh stores mixture thermodynamic properties
    ## 0,  1,     2, 3, 4..4 + ns
    ## cp, gamma, c, e, hi1..hins

    // Compute rho
    fpdtype_t rho = ${" + ".join([f"u[{n}]" for n in range(ns)])};
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${Eix}];

    // Mass fractions
% for n in range(ns):
    q[${n}] = u[${n}]*invrho;
% endfor

    // Compute velocities
% for i in range(ndims):
    q[${i + vix}] = u[${i + vix}]*invrho;
% endfor

    // Compute mixture properties
    fpdtype_t R = 0.0;
% for n in range(ns):
    R += q[${n}]*${c['Ru']/c['MW'][n]};
% endfor

    // Internal energy (per mass)
    fpdtype_t e = rhoE*invrho - 0.5*${pyfr.dot('q[{i}]', i=(vix,vix + ndims))};

    // Iterate on T
    fpdtype_t cp;
    % if fast_props:
    // Quadratic guess for T
    fpdtype_t a = -(${'+'.join([f'{N7[n,1]*Ru/MW[n]/2.0}*q[{n}]' for n in range(ns)])});
    fpdtype_t b = R - (${'+'.join([f'{N7[n,0]*Ru/MW[n]}*q[{n}]' for n in range(ns)])});
    fpdtype_t c = e - (${'+'.join([f'{N7[n,5]*Ru/MW[n]}*q[{n}]' for n in range(ns)])});
    fpdtype_t T = fmin(${c['Tmax']}, fmax(${c['Tmin']}, fabs(a) < ${fpdtype_eps} ? -c/b : (-b + sqrt(fmax(0.0,b*b-4*a*c)))/(2*a)));
    % else:
    fpdtype_t T = ${0.5*(c['Tmax']-c['Tmin'])}; // Initial guess
    % endif
    ${pyfr.expand('T_iter', 'e', 'cp', 'R', 'T', 'q', 'qh')};

    // Equilibrium T, p
    q[${rhoix}] = rho;
    q[${pix}] = rho * R * T;
    q[${Tix}] = T;

    // Mixture gamma, cp
    qh[0] = cp / (cp - R);
    qh[1] = cp;

    // Mixture speed of sound
    qh[2] = sqrt(qh[0] * R * T);

    // internal energy
    qh[3] = rho*e;

    // Store species enthalpy (per mass)
    // ^ done in T_iter

#ifdef DEBUG
  printf("*********************************\n");
  printf("TPG THERMODYNAMIC PROPERTIES\n");
  printf("INPUT STATE\n");
  printf("therm&rho = %e\n", q[${rhoix}]);
  printf("therm&e = %e\n", qh[3]/rho);
% for n in range(ns):
  printf("therm&rhoY_${c['names'][n]} = %e\n", u[${n}]);
% endfor

  printf("\nCOMPUTED STATE\n");
  printf("therm&p = %e\n", q[${pix}]);
  printf("therm&T = %e\n", q[${Tix}]);
% for n in range(ns):
  printf("therm&Y_${c['names'][n]} = %e\n", q[${n}]);
% endfor

  printf("\nCOMPUTED PROPERTIES\n");
  printf("therm&R = %e\n", R);
  printf("therm&gamma = %e\n", qh[0]);
  printf("therm&cp = %e\n", qh[1]);
  printf("therm&c = %e\n", qh[2]);
% for n in range(ns):
  printf("therm&h_${c['names'][n]} = %e\n", qh[${4 + n}]);
% endfor
  printf("*********************************\n");
#endif
</%pyfr:macro>