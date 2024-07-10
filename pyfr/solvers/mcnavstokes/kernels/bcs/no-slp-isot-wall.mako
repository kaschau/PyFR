<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>
<% ns = c['ns'] %>
<% Yix = ndims+2 %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
    fpdtype_t invrho = 1.0/ul[0];

    // Set right primatives
    qr[0] = ql[0];
    // Set zero velocity so we only compute internal energy
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + 1}] = -ql[${i + 1}] + 2*${c[v]};
% endfor
    // Set Temperature
    qr[${ndims+1}] = ${c['T']};

    // Set species
    qr[${nvars}] = 1.0;
% for n in range(ns-1):
    qr[${Yix + n}] = ul[${Yix + n}]*invrho;
    qr[${nvars}] -= qr[${Yix + n}];
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
    fpdtype_t invrho = 1.0/ul[0];

    // Set right primatives
    qr[0] = ql[0];
    // Set zero velocity so we only compute internal energy
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + 1}] = -ql[${i + 1}] + 2*${c[v]};
% endfor
    // Set Temperature
    qr[${ndims+1}] = ${c['T']};

    // Set species
    qr[${nvars}] = 1.0;
% for n in range(ns-1):
    qr[${Yix + n}] = ul[${Yix + n}]*invrho;
    qr[${nvars}] -= qr[${Yix + n}];
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
