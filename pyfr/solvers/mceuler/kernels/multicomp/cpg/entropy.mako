<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    fpdtype_t T = q[${Tix}];
    fpdtype_t rho = q[${rhoix}];
    e = 0.0;
    // Compute mixture entropy
    % for n in range(ns):
    // ${c['names'][n]} Entropy
    {
      <% cvk = c['cp0'][n] - c['Ru']/c['MW'][n] %>
      <% gammak = c['cp0'][n]/cvk %>
      fpdtype_t rhoYk = rho*q[${n}];
      e += rhoYk > 0.0 ? ${cvk}*rhoYk*log(pow(rhoYk, ${1.0 - gammak})*T) : 0.0;
    }
    % endfor
    e = ((T > 0) && (q[0] > 0)) ? e : ${fpdtype_max};
</%pyfr:macro>
