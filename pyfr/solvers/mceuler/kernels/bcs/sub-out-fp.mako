<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
<% ns = c['ns'] %>
<% Yix = ndims+2 %>

    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('mixture_state', 'ul', 'ql', 'qhl')};

    // set right side primatives
    fpdtype_t qr[${nvars+1}];
    // fix pressure
    qr[0] = ${c['p']};

% for i in range(ndims):
    qr[${i+1}] = ql[${i+1}];
% endfor
    qr[${ndims+1}] = ql[${ndims+1}];
% for n in range(ns):
    qr[${Yix+n}] = ql[${Yix+n}];
% endfor

${pyfr.expand('prims_to_cons', 'qr', 'ur')};

</%pyfr:macro>