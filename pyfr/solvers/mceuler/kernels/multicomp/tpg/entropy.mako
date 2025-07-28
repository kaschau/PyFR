<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = 'fast_coeff' in c %>\
% if fast_props:
<% fast_coeff = c['fast_coeff'] %>\
% else:
<% T_cutoff = c['T_cutoff'] %>\
<% NASA7_Thigh = c['NASA7_Thigh'] %>\
<% NASA7_Tlow = c['NASA7_Tlow'] %>\
% endif\

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
        ss = ${pyfr.nasa_s(fast_coeff[n], Ru, MW[n])};
      % else:
        if (T < ${T_cutoff[n]})
        {
          ss = ${pyfr.nasa_s(NASA7_Tlow[n], Ru, MW[n])};
        }else{
          ss = ${pyfr.nasa_s(NASA7_Thigh[n], Ru, MW[n])};
        }
      % endif

      <% Rk = c['Ru']/c['MW'][n] %>\
      s += q[${n}] > 0.0 ? q[${n}] * (ss - ${Rk}*log(u[${n}])) : 0.0;
    }
% endfor

    // Return the specific thermodynamic entropy (mass basis)
    s = (T > 0) ? s : ${-fpdtype_max};
</%pyfr:macro>