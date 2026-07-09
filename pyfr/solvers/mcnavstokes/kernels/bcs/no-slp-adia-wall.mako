<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor
% for i in range(ndims):
    ur[${ns + i}] = -ul[${ns + i}];
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t rho_wla = ${' + '.join(f'ul[{n}]' for n in range(ns))};

% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor
% for i in range(ndims):
    ur[${ns + i}] = 0.0;
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}]
                     - (0.5/rho_wla)*${pyfr.dot('ul[{i}]',
                                                i=(ns, ns + ndims))};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ub, nl, grad_ul, grad_ur'>
    ${fluid.decl('ub', 'rho, invrho, v, Y, T, e_sp', suffix='g')}
    fpdtype_t E_g = ub[${nvars - 1}]*invrhog;

% for d in range(ndims):
    fpdtype_t rho_${'xyz'[d]} = ${' + '.join(f'grad_ul[{d}][{n}]'
                                             for n in range(ns))};
% endfor

    // Velocity derivatives (rho-scaled)
% for i, d in pyfr.ndrange(ndims, ndims):
    fpdtype_t v${i}_${'xyz'[d]} = grad_ul[${d}][${ns + i}]
                                - vg[${i}]*rho_${'xyz'[d]};
% endfor

    // Temperature part of the rhoE gradient (rho*cv*dT/dx scaling)
% for d in range(ndims):
    fpdtype_t eyy_${'xyz'[d]} = ${' + '.join(
        f'e_spg[{n}]*(grad_ul[{d}][{n}] - Yg[{n}]*rho_{"xyz"[d]})'
        for n in range(ns))};
    fpdtype_t T_${'xyz'[d]} = grad_ul[${d}][${nvars - 1}]
        - E_g*rho_${'xyz'[d]}
        - (${' + '.join(f'vg[{i}]*v{i}_{"xyz"[d]}' for i in range(ndims))})
        - eyy_${'xyz'[d]};
% endfor

    // Copy all fluid-side gradients across to wall-side gradients
    ${pyfr.expand('bc_common_grad_copy', 'ub', 'nl', 'grad_ul', 'grad_ur')};

    // Remove the normal component of the temperature gradient
    fpdtype_t Tdotn = ${' + '.join(f'T_{"xyz"[d]}*nl[{d}]'
                                   for d in range(ndims))};
% for d in range(ndims):
    grad_ur[${d}][${nvars - 1}] -= Tdotn*nl[${d}];
% endfor

    // Remove the normal component of each species gradient
    fpdtype_t Ydotn;
% for n in range(ns):
    Ydotn = ${' + '.join(f'(grad_ul[{d}][{n}] - Yg[{n}]*rho_{"xyz"[d]})'
                         f'*nl[{d}]' for d in range(ndims))};
% for d in range(ndims):
    grad_ur[${d}][${n}] -= Ydotn*nl[${d}];
% endfor
% endfor
</%pyfr:macro>
