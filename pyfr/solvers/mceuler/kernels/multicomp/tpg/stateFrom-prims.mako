<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='stateFrom-prims' params='u, q, qh'>
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

    // Compute density
    fpdtype_t rho = q[0]/(R*q[${ndims+1}]);
    u[0] = rho;

    // Compute momentum
% for i in range(ndims):
    u[${i+1}] = q[${i+1}]*rho;
% endfor

    // Total energy
    fpdtype_t T = q[${ndims+1}];
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
    qh[${4+n}] = hs;
    }
% endfor

    fpdtype_t rhoe = rho*h - q[0];
    u[${ndims+1}] = rhoe + 0.5*rho*${pyfr.dot('q[{i}]', i=(1,ndims+1))};

    // Species mass
% for n in range(ns-1):
    u[${Yix+n}] = q[${Yix+n}]*rho;
% endfor

    // Store gamma, cp, c, rhoe
    qh[0] = cp/(cp-R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0]*R*q[${ndims+1}]);
    qh[3] = rhoe

</%pyfr:macro>