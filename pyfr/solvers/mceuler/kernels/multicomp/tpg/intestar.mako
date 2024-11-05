<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = N7.shape[1] == 7 %>\

<%pyfr:macro name='compute_intestar' params='u, q, qh, intestar'>

    fpdtype_t T = q[${Tix}];

    intestar = qh[3]/q[${rhoix}];
    // Compute shifted internal energy
% for n in range(ns):

    % if fast_props:
      intestar -= q[${n}] * ${N7[n,5] * Ru/MW[n]};
    % else:
      if (T < ${N7[n,0]})
      {
        intestar -= q[${n}] * ${N7[n,8 + 5] * Ru/MW[n]};
      }else
      {
        intestar -= q[${n}] * ${N7[n,1 + 5] * Ru/MW[n]};
      }
    % endif
% endfor
</%pyfr:macro>