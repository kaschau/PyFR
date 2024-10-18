<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% div = [1.0, 2.0, 3.0, 4.0, 5.0] %>\

<%pyfr:macro name='stateFrom-prims' params='u, q, qh'>

    ## q is an array of length nvars + 2
    ## storing all primatives
    ## 0:ns-1,    ns:ns+ndims, ns+ndims+1, ns+ndims+2, nvars + 2
    ## Y0...Ynsp, u,v(,w),  rho          , p,          T

    ## qh stores mixture thermodynamic properties
    ## 0,  1,     2, 3, 4..4 + ns
    ## cp, gamma, c, e, hi1..hins


    // Compute mixture properties
    fpdtype_t R = 0.0;
% for n in range(ns):
    R += q[${n}]*${1.0/MW[n]};
% endfor
    R *= ${Ru};

    // Compute density
    fpdtype_t rho = q[${pix}]/(R*q[${Tix}]);
    q[$[rhoix]] = rho;

    // Species mass
% for n in range(ns):
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
    h += hs * q[${n}];
    cp += cps * q[${n}];
    qh[${4 + n}] = hs;
    }
% endfor

    fpdtype_t rhoe = rho*h - q[${rhoix}];
    u[${Eix}] = rhoe + 0.5*rho*${pyfr.dot('q[{i}]', i=(vix,vix + ndims))};

    // Store gamma, cp, c, rhoe
    qh[0] = cp / (cp - R);
    qh[1] = cp;
    qh[2] = sqrt(qh[0]*R*q[${Tix}]);
    qh[3] = rhoe;

    // Store species enthalpy (per mass)
    // ^ done up there

</%pyfr:macro>