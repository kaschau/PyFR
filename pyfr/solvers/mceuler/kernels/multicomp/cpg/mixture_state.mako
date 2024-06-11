<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='mixture_state' params='u, q, qh'>
    fpdtype_t rho = u[0];
    fpdtype_t invrho = 1.0/u[0];
    fpdtype_t rhoE = u[${ndims + 1}];

    q[${ns+1}] = 1.0;
% for n in range(ns-1):
    q[${2+n}] = u[${n+ndims+2}]*invrho;
    q[${ns+1}] -= q[${2+n}];
% endfor

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    R += q[${2+n}]*${1.0/props['MW'][n]};
    cp += q[${2+n}]*${props['cp0'][n]};
% endfor
    R *= ${props['Ru']};

    // Mixture gamma, cp
    qh[0] = cp / (cp - R);
    qh[1] = cp;

    // Internal energy (per mass)
    fpdtype_t e = (rhoE - 0.5*invrho*${pyfr.dot('u[{i}]', i=(1,ndims+1))})*invrho;

    // Equilibrium T, p
    q[1] = e / (cp - R);
    q[0] = rho * R * q[1];

    // Mixture speed of sound
    qh[2] = sqrt(qh[0] * R * q[1]);

</%pyfr:macro>