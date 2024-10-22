<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    e = 0.0;
    // Compute mixture entropy
    % for n in range(ns):
    // ${c['names'][n]} Entropy
    {
      <% Rk = c['Ru']/c['MW'][n] %>\
      <% cvk = c['cp0'][n] - Rk %>\
      <% gammak = c['cp0'][n]/cvk %>\
      fpdtype_t rhoYk = u[${n}];
      fpdtype_t Yk = q[${n}];
      e += (rhoYk > 0.0) ? rhoYk*(${cvk}*log(q[${Tix}]) - ${Rk}*log(rhoYk)) : 0.0;
    }
    % endfor

    // Return the specific thermodynamic entropy (mass basis)
    e = (q[${Tix}] > 0.0) ? e : ${-fpdtype_max};
</%pyfr:macro>
