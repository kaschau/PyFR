<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.mceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.flux'/>

<% tau = c['ldg-tau'] %>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, artviscl, nl, magnl'>

<% ns = c['ns'] %>
    // Compute left thermodynamic quantities
    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('mixture_state', 'ul', 'ql', 'qhl')};

    // Viscous states
    fpdtype_t ur[${nvars}], gradur[${ndims}][${nvars}];
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur')};
    ${pyfr.expand('bc_ldg_grad_state', 'ul', 'nl', 'gradul', 'gradur')};

    // Compute right thermodynamic quantities
    fpdtype_t qr[${nvars+1}];
    fpdtype_t qhr[${3+ns}];
    ${pyfr.expand('mixture_state', 'ur', 'qr', 'qhr')};

    fpdtype_t fvr[${ndims}][${nvars}] = {{0}};
    fpdtype_t qtr[${nvars+1}];
    ${pyfr.expand('mixture_transport', 'ur', 'qr', 'qhr', 'qtr')};
    ${pyfr.expand('viscous_flux_add', 'ur', 'gradur', 'qr', 'qhr', 'qtr', 'fvr')};
    ${pyfr.expand('artificial_viscosity_add', 'gradur', 'fvr', 'artviscl')};

    // Inviscid (Riemann solve) state
    ${pyfr.expand('bc_rsolve_state', 'ul', 'nl', 'ur')};

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
