<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.flux'/>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, nl, magnl'>
    // Viscous states
    fpdtype_t ur[${nvars}], gradur[${ndims}][${nvars}];
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur')};
    ${pyfr.expand('bc_ldg_grad_state', 'ur', 'nl', 'gradul', 'gradur')};

    ${fluid.decl('ur', 'rho, invrho, v, Y, T, cpmix, gammamix, mu, kappa, D, h, e_sp', suffix='vr')}
    fpdtype_t fvr[${ndims}][${nvars}] = {{0}};
    ${pyfr.expand('viscous_flux_add', 'ur', 'gradur', 'rhovr', 'invrhovr',
                  'vvr', 'Yvr', 'Tvr', 'cpmixvr', 'gammamixvr', 'muvr',
                  'kappavr', 'Dvr', 'hvr', 'e_spvr', 'fvr')};

    // Inviscid (Riemann solve) state
    ${pyfr.expand('bc_rsolve_state', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}], fvcomm;
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    fvcomm = ${' + '.join(f'nl[{j}]*fvr[{j}][{i}]' for j in range(ndims))};
    ul[${i}] = magnl*(ficomm[${i}] + fvcomm);
% endfor
</%pyfr:macro>
