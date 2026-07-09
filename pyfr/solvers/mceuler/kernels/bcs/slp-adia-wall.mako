<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t nor = ${' + '.join(f'ul[{ns + i}]*nl[{i}]'
                                 for i in range(ndims))};

% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor
% for i in range(ndims):
    ur[${ns + i}] = ul[${ns + i}] - 2*nor*nl[${i}];
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>
