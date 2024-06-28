<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t nor = ${' + '.join(f'ul[{i + 1}]*nl[{i}]' for i in range(ndims))};
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor
    ur[${ndims + 1}] = ul[${ndims + 1}];
% for n in range(c['ns']-1):
    ur[${ndims + 2 + n}] = ul[${ndims + 2 + n}];
% endfor
</%pyfr:macro>
