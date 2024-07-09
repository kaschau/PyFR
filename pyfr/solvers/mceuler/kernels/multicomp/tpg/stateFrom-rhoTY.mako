<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## Here we assume q is complete, except for pressure
## Assume u is incomplete, except density

<%pyfr:macro name='stateFrom-rhoTY' params='u, q, qh'>
<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\
<% div = [1.0, 2.0, 3.0, 4.0, 5.0] %>\

    ## q is an array of length nvars + 1
    ## storing all primatives
    ## 0, 1:ndims, ndims+1, ndims+2 ... nvars
    ## p, u,v(,w),   T    ,   Y0    ...  Yns

    // Compute mixture properties
    fpdtype_t R = 0.0;
% for n in range(ns):
    R += q[${Yix+n}]*${1.0/MW[n]};
% endfor
    R *= ${Ru};

    fpdtype_t T = q[${ndims+1}];
    fpdtype_t rho = u[0];

    // Pressure
    q[0] = rho*R*T;

    // Compute momentum
% for i in range(ndims):
    u[${i+1}] = q[${i+1}]*rho;
% endfor

    // Total Energy
    fpdtype_t h = 0.0;
    fpdtype_t cp = 0.0;
% for n in range(ns):
    // ${c['names'][n]} Properties
    {
    fpdtype_t hs, cps;
    if (T < ${N7[n,0]})
    {
    <% m = 8 %>
        cps = ${'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n])+')'*4};
        hs = T*(${'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n]/div)+')'*4}) + ${N7[n, m + 5] * Ru/MW[n]};
    }else
    {
    <% m = 1 %>
        cps = ${'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n])+')'*4};
        hs = T*(${'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n]/div)+')'*4}) + ${N7[n, m + 5] * Ru/MW[n]};
    }
    h += hs * q[${Yix+n}];
    cp += cps * q[${Yix+n}];
    qh[${3+n}] = hs;
    }
% endfor

    // Total Energy
    fpdtype_t e = h - R*T;
    u[${ndims+1}] = rho*(e + 0.5*${pyfr.dot('q[{i}]', i=(1,ndims+1))});

% for n in range(ns-1):
    u[${Yix+n}] = rho*q[${Yix+n}];
% endfor

    // Store gamma, cp, c
    qh[0] = cp/(cp-R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0]*R*T);

</%pyfr:macro>