<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## Here we assume q is complete, except for pressure
## Assume u is incomplete, except density

<%pyfr:macro name='stateFrom-rhoTY' params='u, q, qh'>
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    R += q[${Yix+n}]*${1.0/c['MW'][n]};
    cp += q[${Yix+n}]*${c['cp0'][n]};
% endfor
    R *= ${c['Ru']};

    fpdtype_t T = q[${ndims+1}];
    fpdtype_t rho = u[0];

    // Pressure
    q[0] = rho*R*T;

    // Compute momentum
% for i in range(ndims):
    u[${i+1}] = q[${i+1}]*rho;
% endfor

    // Total Energy
    fpdtype_t e = (cp-R) * T;
    u[${ndims+1}] = rho*(e + 0.5*${pyfr.dot('q[{i}]', i=(1,ndims+1))});

% for n in range(ns-1):
    u[${Yix+n}] = rho*q[${Yix+n}];
% endfor

    // Store gamma, cp, c, hs
    qh[0] = cp/(cp-R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0]*R*T);
% for n in range(ns):
    qh[${3+n}] = T*${c['cp0'][n]};
% endfor

</%pyfr:macro>