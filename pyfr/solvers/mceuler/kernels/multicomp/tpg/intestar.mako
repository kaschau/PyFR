<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\

<%pyfr:macro name='compute_intestar' params='u, q, qh, intestar'>

    fpdtype_t T = q[${Tix}];

    intestar = qh[3]/q[${rhoix}];
    // Compute shifted internal energy
% for n in range(ns):
      if (T < ${N7[n,0]})
      {
      <% m = 8 %>
          intestar -= q[${n}] * ${N7[n,m + 5] * Ru/MW[n]};
      }else
      {
      <% m = 1 %>
          intestar -= q[${n}] * ${N7[n,m + 5] * Ru/MW[n]};
      }
% endfor
</%pyfr:macro>