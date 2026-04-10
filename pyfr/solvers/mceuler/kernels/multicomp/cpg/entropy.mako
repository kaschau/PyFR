<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='compute_entropy' params='u, q, s'>

    s = 0.0;
    % for n in range(mcf.ns):
    {
      <% Rk = mcf.Ru / mcf[n].MW %>
      <% cvk = mcf[n].cp0 - Rk %>
      fpdtype_t Yk = q[${n}];
      s += (Yk > 0.0) ? Yk*(${cvk}*log(q[${Tix}]) - ${Rk}*log(u[${n}])) : 0.0;
    }
    % endfor

    s = (q[${Tix}] > 0.0) ? s : ${-fpdtype_max};
</%pyfr:macro>
