<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    fpdtype_t T = q[${ndims + 1}];
    e = 0.0;
    // Compute mixture entropy
    % for n in range(ns):
    {
      fpdtype_t cvk = ${c['cp0'][n]} - ${c['Ru']/c['MW'][n]};
      fpdtype_t gammak = ${c['cp0'][n]}/cvk;
      e += cvk*u[${Yix + n}] * log(pow(fmax(1e-6, u[${Yix+n}]), 1.0 - gammak)*T);
    }
    % endfor
    e = exp(e);
</%pyfr:macro>
