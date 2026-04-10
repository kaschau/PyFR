<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='compute_entropy' params='u, q, s'>

    fpdtype_t T = q[${Tix}];
    fpdtype_t logT = log(T);
    s = 0.0;
% for n in range(mcf.ns):
    {
      fpdtype_t ss = ${mcf[n].s_expr('T')};
      <% Rk = mcf.Ru / mcf[n].MW %>
      s += q[${n}] > 0.0 ? q[${n}] * (ss - ${Rk}*log(u[${n}])) : 0.0;
    }
% endfor

    s = (T > 0) ? s : ${-fpdtype_max};
</%pyfr:macro>
