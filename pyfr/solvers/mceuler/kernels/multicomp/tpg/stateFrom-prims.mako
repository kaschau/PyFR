<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='stateFrom-prims' params='u, q, qh'>

    // Compute mixture gas constant
    fpdtype_t R = 0.0;
% for n in range(mcf.ns):
    R += q[${n}]*${mcf.Ru / mcf[n].MW};
% endfor

    // Compute density
    fpdtype_t rho = q[${pix}]/(R*q[${Tix}]);
    q[${rhoix}] = rho;

    // Species mass
% for n in range(mcf.ns):
    u[${n}] = q[${n}]*rho;
% endfor

    // Compute momentum
% for i in range(ndims):
    u[${i + vix}] = q[${i + vix}]*rho;
% endfor

    // Total energy
    fpdtype_t T = q[${Tix}];
    fpdtype_t h = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(mcf.ns):
    {
      fpdtype_t cps = ${mcf[n].cp_expr('T')};
      fpdtype_t hs = ${mcf[n].h_expr('T')};
      h += hs * q[${n}];
      cp += cps * q[${n}];
      qh[${4 + n}] = hs;
    }
% endfor

    fpdtype_t rhoe = rho*h - q[${pix}];
    u[${Eix}] = rhoe + 0.5*rho*${pyfr.dot('q[{i}]', i=(vix,vix + ndims))};

    // Store gamma, cp, c, rhoe
    qh[0] = cp / (cp - R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0]*R*q[${Tix}]);
    qh[3] = rhoe;

</%pyfr:macro>
