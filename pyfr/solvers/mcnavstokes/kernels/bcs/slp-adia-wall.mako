<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.flux'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr'>
    fpdtype_t nor = ${' + '.join(f'ul[{i + vix}]*nl[{i}]' for i in range(ndims))};

    // Species
% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor

    // Momentum
% for i in range(ndims):
    ur[${i + vix}] = ul[${i + vix}] - 2*nor*nl[${i}];
% endfor

    // Total energy
    ur[${Eix}] = ul[${Eix}];
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state' params='ul, ql, qhl, gradul, artviscl, nl, magnl'>

## bc_ldg_state AND bc_rsolve_state must fill in ur, qr, qhr

    // Ghost state r
    fpdtype_t ur[${nvars}];
    fpdtype_t qr[${nvars + 2}];
    fpdtype_t qhr[${4 + ns}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'ql', 'qhl', 'nl', 'ur', 'qr', 'qhr')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'ql', 'qr', 'qhl', 'qhr', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*(ficomm[${i}]);
% endfor
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
