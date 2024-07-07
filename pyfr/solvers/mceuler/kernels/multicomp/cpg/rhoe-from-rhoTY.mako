<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='rhoe-from-rhoTY' params='q, u'>
    ## q is an array of length nvars + 1
    ## storing all primatives
    ## 0, 1:ndims, ndims+1, ndims+2 ... nvars
    ## p, u,v(,w),   T    ,   Y0    ...  Yns

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

    // Internal energy
    u[${ndims+1}] = (cp-R) * q[${ndims+1}]*u[0];

</%pyfr:macro>