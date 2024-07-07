<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='rhoe-from-rhoTY' params='q, u'>
    ## q is an array of length nvars + 1
    ## storing all primatives
    ## 0, 1:ndims, ndims+1, ndims+2 ... nvars
    ## p, u,v(,w),   T    ,   Y0    ...  Yns

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\
<% div = [1.0, 2.0, 3.0, 4.0, 5.0] %>\

    // Compute mixture properties
    fpdtype_t R = 0.0;
% for n in range(ns):
    R += q[${Yix+n}]*${1.0/MW[n]};
% endfor
    R *= ${Ru};

    // Internal energy
    fpdtype_t T = q[${ndims+1}];
    fpdtype_t h = 0.0;
% for n in range(ns):
    // ${c['names'][n]} Properties
    {
    fpdtype_t hs;
    if (T < ${N7[n,0]})
    {
<% m = 8 %>
        hs = T*(${'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n]/div)+')'*4}) + ${N7[n, m + 5] * Ru/MW[n]};
    }else
    {
<% m = 1 %>
        hs = T*(${'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n]/div)+')'*4}) + ${N7[n, m + 5] * Ru/MW[n]};
    }
    h += hs * q[${Yix+n}];
    }

    fpdtype_t rhoe = rho*h - q[0];
    u[${ndims+1}] = rhoe;

</%pyfr:macro>