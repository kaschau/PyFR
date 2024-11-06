<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = N7.shape[1] == 7 %>\

<%pyfr:macro name='T_iter' params='e, cp, Rmix, T, q, qh'>

    <% tol = 1e-8 %>\
    <% niter_max = 1e-8 %>\
    fpdtype_t error = ${fpdtype_max};
    for (int niter = 0; niter < ${niter_max} && abs(error) > ${tol}; niter++)
    {
        fpdtype_t h = 0.0;
        cp = 0.0;
% for n in range(ns):
        // ${c['names'][n]} Properties
        {
        fpdtype_t cps, hs;
        % if fast_props:
            cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 0)};
            hs = ${pyfr.nasa_hs(N7[n,:], Ru, MW[n], 0)};
        % else:
          if (T < ${N7[n,0]})
          {
            cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 8)};
            hs = ${pyfr.nasa_hs(N7[n,:], Ru, MW[n], 8)};
          }else
          {
            cps = ${pyfr.nasa_cps(N7[n,:], Ru, MW[n], 1)};
            hs = ${pyfr.nasa_hs(N7[n,:], Ru, MW[n], 1)};
          }
        % endif
        cp += cps * q[${n}];
        h += hs * q[${n}];
        qh[${4 + n}] = hs;
        }
% endfor
    error = e - (h - Rmix * T);
    // Newton's Method
    T = T - error / (-cp - Rmix);
    }
</%pyfr:macro>