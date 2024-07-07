<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
<% ns = c['ns'] %>
    // Set right thermodynamic quantities
    fpdtype_t qr[${nvars+1}];

    qr[${ndims+1}] = ${c['T']};
    qr[${nvars}] = 1.0;
% for n in range(ns-1):
    qr[${ndims + 2 + n}] = ul[${ndims + 2 + n}]/ul[0];
    qr[${nvars}] -= qr[${ndims + 2 + n}];
% endfor

    ur[0] = ul[0];
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = -ul[${i + 1}] + 2*${c[v]}*ul[0];
% endfor

    ${pyfr.expand('rhoe-from-rhoTY', 'qr', 'ur')};
    ur[${ndims + 1}] += 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};

% for n in range(ns-1):
    ur[${ndims + 2 + n}] = ul[${ndims + 2 + n}];
% endfor

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
<% ns = c['ns'] %>
    // Set right thermodynamic quantities
    fpdtype_t qr[${nvars+1}];

    qr[${ndims+1}] = ${c['T']};
    qr[${nvars}] = 1.0;
% for n in range(ns-1):
    qr[${ndims + 2 + n}] = ul[${ndims + 2 + n}]/ul[0];
    qr[${nvars}] -= qr[${ndims + 2 + n}];
% endfor

    ur[0] = ul[0];
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + 1}] = -ul[${i + 1}] + 2*${c[v]}*ul[0];
% endfor

    ${pyfr.expand('rhoe-from-rhoTY', 'qr', 'ur')};
    ur[${ndims + 1}] += 0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};

% for n in range(ns-1):
    ur[${ndims + 2 + n}] = ul[${ndims + 2 + n}];
% endfor

</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
