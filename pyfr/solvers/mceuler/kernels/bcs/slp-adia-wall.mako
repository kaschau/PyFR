<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    fpdtype_t nor = ${' + '.join(f'ul[{vix + i}]*nl[{i}]' for i in range(ndims))};

% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor

% for i in range(ndims):
    ur[${i + vix}] = ul[${i + vix}] - 2*nor*nl[${i}];
% endfor

    ur[${Eix}] = ul[${Eix}];

    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};

</%pyfr:macro>
