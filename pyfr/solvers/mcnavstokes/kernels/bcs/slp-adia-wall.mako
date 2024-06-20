<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.mixture_state'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.flux'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur'>
    fpdtype_t nor = ${' + '.join(f'ul[{i + 1}]*nl[{i}]' for i in range(ndims))};
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor
    ur[${ndims + 1}] = ul[${ndims + 1}];
% for n in range(ns-1):
    ur[${ndims + 2 + n}] = ul[${ndims + 2 + n}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, artviscl, nl, magnl'>
    // Ghost state r
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'nl', 'ur')};

<% ns = c['ns'] %>
    // Compute left thermodynamic quantities
    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('mixture_state', 'ul', 'ql', 'qhl')};

    // Compute right thermodynamic quantities
    fpdtype_t qr[${nvars+1}];
    fpdtype_t qhr[${3+ns}];
    ${pyfr.expand('mixture_state', 'ur', 'qr', 'qhr')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'ql', 'qr', 'qhl', 'qhr', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*(ficomm[${i}]);
% endfor
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
