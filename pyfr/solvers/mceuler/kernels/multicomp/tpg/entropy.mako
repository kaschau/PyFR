<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% strict = N7.shape[1] == 15 %>\

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    fpdtype_t T = q[${Tix}];
    fpdtype_t logT = log(T);
    e = 0.0;
    // Compute mixture entropy
% for n in range(ns):
    // ${c['names'][n]} Entropy
    {
      fpdtype_t es;
      % if not strict:
        es = ${pyfr.nasa_scs(N7[n,:], Ru, MW[n], 0)};
      % else:
        if (T < ${N7[n,0]})
        {
          es = ${pyfr.nasa_scs(N7[n,:], Ru, MW[n], 8)};
        }else{
          es = ${pyfr.nasa_scs(N7[n,:], Ru, MW[n], 1)};
        }
      % endif

      <% Rk = c['Ru']/c['MW'][n] %>\
      e += u[${n}] > 0.0 ? u[${n}] * (es - Rk*log(u[${n}])) : 0.0;
    }

    // Return the specific thermodynamic entropy (mass basis)
    e = (T > 0) ? e : ${-fpdtype_max};
% endfor
</%pyfr:macro>