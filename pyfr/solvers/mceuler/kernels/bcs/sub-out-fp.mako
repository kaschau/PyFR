<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
<% ns = c['ns'] %>
<% Yix = ndims+2 %>


    // set right side primatives
    // fix pressure
    qr[0] = ${c['p']};

% for i in range(ndims):
    qr[${i+1}] = ql[${i+1}];
% endfor
    qr[${ndims+1}] = ql[${ndims+1}];
% for n in range(ns):
    qr[${Yix+n}] = ql[${Yix+n}];
% endfor

${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>