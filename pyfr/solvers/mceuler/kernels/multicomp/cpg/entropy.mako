<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    fpdtype_t T = q[${ndims + 1}];
    fpdtype_t rho = u[0];
    e = 0.0;
    fpdtype_t Ymin = ${fpdtype_max};
    fpdtype_t Ysum = 0.0;
    // Compute mixture entropy
    % for n in range(ns):
    // ${c['names'][n]} Entropy
    {
      <% cvk = c['cp0'][n] - c['Ru']/c['MW'][n] %>\
      <% gammak = c['cp0'][n]/cvk %>\
      fpdtype_t rhoYk = rho*q[${Yix + n}];
      e += ${cvk}*rhoYk * log(pow(fmax(1e-12, rhoYk), ${1.0 - gammak})*T);
      Ymin = fmin(Ymin, q[${Yix + n}]);
      Ysum += q[${Yix + n}];
    }
    % endfor
    e /= rho;
    e = ((T > 0.0) && (q[0] > 0.0) && (Ymin > 0.0) && (Ysum < 1.0)) ? exp(e) : ${fpdtype_max};
</%pyfr:macro>
