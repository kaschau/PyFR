<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% div = [1.0, 2.0, 3.0, 4.0, 5.0] %>\

<%pyfr:macro name='T_iter' params='e, cp, Rmix, T, q, qh'>

    fpdtype_t tol = 1e-8;
    fpdtype_t error = ${fpdtype_max};
    fpdtype_t niter = 0;
    while ((abs(error) > tol) && (niter < 50))
    {
        fpdtype_t h = 0.0;
        cp = 0.0;
% for n in range(ns):
        // ${c['names'][n]} Properties
        {
        fpdtype_t cps, hs;
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
        cp += cps * q[${n}];
        h += hs * q[${ix+n}];
        qh[${4 + n}] = hs;
        }
% endfor
    error = e - (h - Rmix * T);
    // Newton's Method
    T = T - error / (-cp - Rmix);
    niter += 1;
    }
</%pyfr:macro>