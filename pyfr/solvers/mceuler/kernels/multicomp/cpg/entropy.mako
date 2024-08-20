<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    fpdtype_t T = q[${ndims + 1}];
    fpdtype_t rho = u[0];
    e = 0.0;
    fpdtype_t rhoYmin = ${fpdtype_max};
    // Compute mixture entropy
    % for n in range(ns):
    // ${c['names'][n]} Entropy
    {
      <% cvk = c['cp0'][n] - c['Ru']/c['MW'][n] %>
      <% gammak = c['cp0'][n]/cvk %>
      fpdtype_t rhoYk = rho*q[${Yix + n}];
      e += rhoYk > 0.0 ? ${cvk}*rhoYk*log(pow(rhoYk, ${1.0 - gammak})*T) : 0.0;
    }
    % endfor
    e = ((T > 0) && (q[0] > 0)) ? e : ${fpdtype_max};
</%pyfr:macro>
