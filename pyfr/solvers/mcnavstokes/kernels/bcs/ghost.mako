<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.mceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.flux'/>

## bc_ldg_state AND bc_rsolve_state must fill in ur, qr, qhr

<% tau = c['ldg-tau'] %>
<% ns = c['ns'] %>

<%pyfr:macro name='bc_common_flux_state' params='ul, ql, qhl, gradul, artviscl, nl, magnl'>

    // Right states
    fpdtype_t ur[${nvars}], gradur[${ndims}][${nvars}];
    fpdtype_t qr[${nvars+1}];
    fpdtype_t qhr[${3+ns}];

    ${pyfr.expand('bc_ldg_state', 'ul', 'ql', 'qhl', 'nl', 'ur', 'qr', 'qhr')};
    ${pyfr.expand('bc_ldg_grad_state', 'ul', 'ql', 'qhl', 'nl', 'gradul', 'gradur')};

    // Mixture transport properties
    fpdtype_t qtr[${nvars+1}];
    ${pyfr.expand('mixture_transport', 'ur', 'qr', 'qhr', 'qtr')};

    fpdtype_t fvr[${ndims}][${nvars}] = {{0}};
    ${pyfr.expand('viscous_flux_add', 'ur', 'gradur', 'qr', 'qhr', 'qtr', 'fvr')};
    ${pyfr.expand('artificial_viscosity_add', 'gradur', 'fvr', 'artviscl')};

    // Inviscid (Riemann solve) state
    ${pyfr.expand('bc_rsolve_state', 'ul', 'ql', 'qhl', 'nl', 'ur', 'qr', 'qhr')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}], fvcomm;
    ${pyfr.expand('rsolve', 'ul', 'ur', 'ql', 'qr', 'qhl', 'qhr', 'nl', 'ficomm')};

% for i in range(nvars):
    fvcomm = ${' + '.join(f'nl[{j}]*fvr[{j}][{i}]' for j in range(ndims))};
% if tau != 0.0:
    fvcomm += ${tau}*(ul[${i}] - ur[${i}]);
% endif

    ul[${i}] = magnl*(ficomm[${i}] + fvcomm);
% endfor
</%pyfr:macro>
