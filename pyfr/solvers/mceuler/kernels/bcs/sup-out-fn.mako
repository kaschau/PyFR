<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
% for i in range(nvars):
    ur[${i}] = ul[${i}];
% endfor
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>
