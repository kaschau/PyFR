<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.mixture_state'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.prims_to_cons'/>

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
% endif
    qr[${ndims+1}] = ql[${ndims+1}];
% for n in range(ns):
    qr[${Yix+n}] = ql[${Yix+n}];
% endfor

${pyfr.expand('prims_to_cons', 'qr', 'ur')};

</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
