<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='T_iter' params='e, h, cp, Rmix, T, q, qh'>
<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\
<% div = [1.0, 2.0, 3.0, 4.0, 5.0] %>\

    fpdtype_t tol = 1e-8;
    fpdtype_t error = 1e100;
    while (abs(error) > tol)
    {
        h = 0.0;
        cp = 0.0;
% for n in range(ns):
        // ${c['names'][n]} Properties
        if (T < ${N7[n,0]})
        {
<% m = 8 %>
            fpdtype_t cps = ${f'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n])+')'*4};
            qh[${3+n}] = T*(${f'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n]/div)+')'*4}) + ${N7[n, m + 5] * Ru/MW[n]};
            cp += cps * q[${Yix+n}];
            h += qh[${3+n}] * q[${Yix+n}];
        }else
        {
<% m = 1 %>
            fpdtype_t cps = ${f'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n])+')'*4};
            qh[${3+n}] = T*(${f'+ T*('.join(str(c) for c in N7[n,m:m+5]*Ru/MW[n]/div)+')'*4}) + ${N7[n, m + 5] * Ru/MW[n]};
            cp += cps * q[${Yix+n}];
            h += qh[${3+n}] * q[${Yix+n}];
        }
% endfor
    error = e - (h - Rmix * T);
    // Newton's Method
    T = T - error / (-cp - Rmix);
    }
</%pyfr:macro>