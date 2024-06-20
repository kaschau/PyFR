<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='T_iter' params='e, h, cp, Rmix, T, q, qh'>

<% N7 = c['NASA7'] %>
<% Ru = c['Ru'] %>
<% MW = c['MW'] %>
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>

    for (ixdtype_t i = 0; i < 50; i++)
    {
        h = 0.0;
        cp = 0.0;
        fpdtype_t T2 = T*T;
        fpdtype_t T3 = T2*T;
        fpdtype_t T4 = T3*T;
        fpdtype_t T5 = T4*T;
        fpdtype_t T2o2 = T2 / 2.0;
        fpdtype_t T3o3 = T3 / 3.0;
        fpdtype_t T4o4 = T4 / 4.0;
        fpdtype_t T5o5 = T5 / 5.0;

        ixdtype_t m;
        fpdtype_t cps;

        if (T < ${N7[0,0]})
        {
<% m = 8 %>
% for n in range(ns):
            // ${c['names'][n]}
            cps = (${N7[n, m + 0] * Ru/MW[n]} +
                   ${N7[n, m + 1] * Ru/MW[n]} * T +
                   ${N7[n, m + 2] * Ru/MW[n]} * T2 +
                   ${N7[n, m + 3] * Ru/MW[n]} * T3 +
                   ${N7[n, m + 4] * Ru/MW[n]} * T4);

            qh[${3+n}] = (${N7[n, m + 0] * Ru/MW[n]} * T +
                          ${N7[n, m + 1] * Ru/MW[n]} * T2o2 +
                          ${N7[n, m + 2] * Ru/MW[n]} * T3o3 +
                          ${N7[n, m + 3] * Ru/MW[n]} * T4o4 +
                          ${N7[n, m + 4] * Ru/MW[n]} * T5o5 +
                          ${N7[n, m + 5] * Ru/MW[n]});
            cp += cps * q[${Yix+n}];
            h += qh[${3+n}] * q[${Yix+n}];
% endfor
        }else
        {
<% m = 1 %>
% for n in range(ns):
            // ${c['names'][n]}
            cps = (${N7[n, m + 0]} +
                   ${N7[n, m + 1] * Ru/MW[n]} * T +
                   ${N7[n, m + 2] * Ru/MW[n]} * T2 +
                   ${N7[n, m + 3] * Ru/MW[n]} * T3 +
                   ${N7[n, m + 4] * Ru/MW[n]} * T4);

            qh[${3+n}] = (${N7[n, m + 0] * Ru/MW[n]} * T +
                          ${N7[n, m + 1] * Ru/MW[n]} * T2o2 +
                          ${N7[n, m + 2] * Ru/MW[n]} * T3o3 +
                          ${N7[n, m + 3] * Ru/MW[n]} * T4o4 +
                          ${N7[n, m + 4] * Ru/MW[n]} * T5o5 +
                          ${N7[n, m + 5] * Ru/MW[n]});
            cp += cps * q[${Yix+n}];
            h += qh[${3+n}] * q[${Yix+n}];
% endfor
        }
    fpdtype_t error = e - (h - Rmix * T);
    T = T - error / (-cp - Rmix);
    }
</%pyfr:macro>