<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\

    fpdtype_t nor = ${' + '.join(f'ul[{i + 1}]*nl[{i}]' for i in range(ndims))};

    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor
    ur[${ndims + 1}] = ul[${ndims + 1}];
% for n in range(ns-1):
    ur[${Yix + n}] = ul[${Yix + n}];
% endfor

    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};

</%pyfr:macro>
