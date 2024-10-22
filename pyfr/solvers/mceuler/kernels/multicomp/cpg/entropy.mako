<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    e = 0.0;
    // Compute mixture entropy
    % for n in range(ns):
    // ${c['names'][n]} Entropy
    {
      <% cvk = c['cp0'][n] - c['Ru']/c['MW'][n] %>
      <% gammak = c['cp0'][n]/cvk %>
      fpdtype_t rhoYk = u[${n}];
      fpdtype_t Yk = q[${n}];
      e += rhoYk > 0.0 ? ${cvk}*Yk*log(pow(rhoYk, ${1.0 - gammak})*q[${Tix}]) : 0.0;
    }
    % endfor

    // Return the specific thermodynamic entropy (mass basis)
    e = (q[${Tix}] > 0) ? e : ${fpdtype_max};
</%pyfr:macro>
