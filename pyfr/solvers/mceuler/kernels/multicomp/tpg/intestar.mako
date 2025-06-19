<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% fast_props = 'fast_coeff' in c %>\
% if fast_props:
<% fast_coeff = c['fast_coeff'] %>\
% else:
<% T_cutoff = c['T_cutoff'] %>\
<% NASA7_Thigh = c['NASA7_Thigh'] %>\
<% NASA7_Tlow = c['NASA7_Tlow'] %>\
% endif\

<%pyfr:macro name='compute_intestar' params='u, q, qh, intestar'>

    fpdtype_t T = q[${Tix}];

    intestar = qh[3];
    // Compute shifted internal energy
% for n in range(ns):

    % if fast_props:
      intestar -= u[${n}] * ${fast_coeff[n][-2] * Ru/MW[n]};
    % else:
      if (T < ${T_cutoff[n]})
      {
        intestar -= u[${n}] * ${NASA7_Tlow[n][5] * Ru/MW[n]};
      }else
      {
        intestar -= u[${n}] * ${NASA7_Thigh[n][5] * Ru/MW[n]};
      }
    % endif
% endfor
</%pyfr:macro>