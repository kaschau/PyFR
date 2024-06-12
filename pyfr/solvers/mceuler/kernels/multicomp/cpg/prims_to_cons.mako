<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='prims_to_cons' params='q, u'>
    ## q is an array of length nvars + 1
    ## storing all primatives
    ## 0, 1:ndims, ndims+1, ndims+2 ... nvars
    ## p, u,v(,w),   T    ,   Y0    ...  Yns

<% Yix = ndims + 2 %>

    // Compute mixture properties
    fpdtype_t R = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    R += q[${Yix+n}]*${1.0/props['MW'][n]};
    cp += q[${Yix+n}]*${props['cp0'][n]};
% endfor
    R *= ${props['Ru']};

    // Compute density
    fpdtype_t rho = q[0]/(R*q[{ndims+1}]);

    // Compute momentum
% for i in range(ndims):
    u[${i+1}] = q[${i+1}]*rho;
% endfor

    // Total energy
    fpdtype_t e = cp/(cp-R) * T;
    u[${ndims+1}] = rho*e + 0.5*rho*${pyfr.dot('q[{i}]', i=(1,ndims+1))};

    // Species mass
% for n in range(ns-1):
    u[${Yix+n}] = q[${Yix+n}]*rho;
% endfor

</%pyfr:macro>