<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
<% ns = c['ns'] %>
    // Compute left thermodynamic quantities
    //HERE
    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('mixture_state', 'ul', 'ql', 'qhl')};

    // Set right thermodynamic quantities
    fpdtype_t qr[${nvars+1}];

    qr[0] = ql[0];
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + 1}] = -ql[${i + 1}] + 2*${c[v]};
% endfor
    qr[${ndims+1}] = ${c['T']};
% for n in range(ns):
    qr[${ndims + 2 + n}] = ql[${ndims + 2 + n}];
% endfor

    ${pyfr.expand('prims_to_cons', 'qr', 'ur')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
<% ns = c['ns'] %>
    // Compute left thermodynamic quantities
    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('mixture_state', 'ul', 'ql', 'qhl')};

    // Set right thermodynamic quantities
    fpdtype_t qr[${nvars+1}];

    qr[0] = ql[0];
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + 1}] = ql[${i + 1}];
% endfor
    qr[${ndims+1}] = ${c['T']};
% for n in range(ns):
    qr[${ndims + 2 + n}] = ql[${ndims + 2 + n}];
% endfor

    ${pyfr.expand('prims_to_cons', 'qr', 'ur')};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
