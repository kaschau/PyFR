<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='mixture_state' params='u, q, qh'>
    ## q is an array of length nvars + 1
    ## storing all primatives
    ## 0, 1:ndims, ndims+1, ndims+2 ... nvars
    ## p, u,v(,w),   T    ,   Y0    ...  Yns
    fpdtype_t rho = u[0];
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${ndims + 1}];

    // Compute velocities
% for i in range(ndims):
    q[${i+1}] = u[${i+1}]*invrho;
% endfor

## First species index
<% Yix = ndims + 2 %>

    q[${nvars}] = 1.0;
% for n in range(ns-1):
    q[${Yix+n}] = u[${Yix+n}]*invrho;
    q[${nvars}] -= q[${Yix+n}];
% endfor

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    R += q[${Yix+n}]*${1.0/props['MW'][n]};
    cp += q[${Yix+n}]*${props['cp0'][n]};
% endfor
    R *= ${props['Ru']};

    // Mixture gamma, cp
    qh[0] = cp / (cp - R);
    qh[1] = cp;

    // Internal energy (per mass)
    fpdtype_t e = (rhoE - 0.5*rho*${pyfr.dot('q[{i}]', i=(1,ndims+1))})*invrho;

    // Equilibrium T, p
    q[${ndims + 1}] = e / (cp - R);
    q[0] = rho * R * q[${ndims + 1}];

    // Mixture speed of sound
    qh[2] = sqrt(qh[0] * R * q[${ndims + 1}]);

</%pyfr:macro>