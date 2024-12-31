<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = N7.shape[1] == 7 %>\

<%pyfr:macro name='compute_entropy' params='u, q, s'>

    fpdtype_t T = q[${Tix}];
    fpdtype_t logT = log(T);
    s = 0.0;
    // Compute mixture entropy
% for n in range(ns):
    // ${c['names'][n]} Entropy
    {
      fpdtype_t ss;
      % if fast_props:
        ss = ${pyfr.nasa_s(N7[n,:], Ru, MW[n], 0)};
      % else:
        if (T < ${N7[n,0]})
        {
          ss = ${pyfr.nasa_s(N7[n,:], Ru, MW[n], 8)};
        }else{
          ss = ${pyfr.nasa_s(N7[n,:], Ru, MW[n], 1)};
        }
      % endif

      s += q[${n}] > 0.0 ? q[${n}] * ss : 0.0;
    }

    // Return the specific thermodynamic entropy (mass basis)
    s = (T > 0) ? s : ${-fpdtype_max};
% endfor
</%pyfr:macro>