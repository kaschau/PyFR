<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ${fluid.decl('ul', 'rho, invrho, p, Y, v', suffix='wr')}

    fpdtype_t T_wr = ${c['T']};
    fpdtype_t R_wr = ${' + '.join(f'Ywr[{n}]*{r}'
                                  for n, r in enumerate(sp_R))};
    fpdtype_t h_wr = ${' + '.join(f"Ywr[{n}]*({fluid.h_of_T(n, 'T_wr')})"
                                  for n in range(ns))};
    fpdtype_t rho_wr = pwr/(R_wr*T_wr);
    fpdtype_t vb_wr[${ndims}];
% for i, v in enumerate('uvw'[:ndims]):
    vb_wr[${i}] = -vwr[${i}] + 2*(${c[v]});
% endfor

% for n in range(ns):
    ur[${n}] = rho_wr*Ywr[${n}];
% endfor
% for i in range(ndims):
    ur[${ns + i}] = rho_wr*vb_wr[${i}];
% endfor
    ur[${nvars - 1}] = rho_wr*h_wr - pwr
                     + 0.5*rho_wr*${pyfr.dot('vb_wr[{i}]', i=ndims)};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
    ${fluid.decl('ul', 'rho, invrho, p, Y, v', suffix='wl')}

    fpdtype_t T_wl = ${c['T']};
    fpdtype_t R_wl = ${' + '.join(f'Ywl[{n}]*{r}'
                                  for n, r in enumerate(sp_R))};
    fpdtype_t h_wl = ${' + '.join(f"Ywl[{n}]*({fluid.h_of_T(n, 'T_wl')})"
                                  for n in range(ns))};
    fpdtype_t rho_wl = pwl/(R_wl*T_wl);
    fpdtype_t vb_wl[${ndims}];
% for i, v in enumerate('uvw'[:ndims]):
    vb_wl[${i}] = -vwl[${i}] + 2*(${c[v]});
% endfor

% for n in range(ns):
    ur[${n}] = rho_wl*Ywl[${n}];
% endfor
% for i in range(ndims):
    ur[${ns + i}] = rho_wl*vb_wl[${i}];
% endfor
    ur[${nvars - 1}] = rho_wl*h_wl - pwl
                     + 0.5*rho_wl*${pyfr.dot('vb_wl[{i}]', i=ndims)};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ub, nl, grad_ul, grad_ur'>
    ${fluid.decl('ub', 'rho, invrho, Y', suffix='gi')}

% for d in range(ndims):
    fpdtype_t rho_${'xyz'[d]} = ${' + '.join(f'grad_ul[{d}][{n}]'
                                             for n in range(ns))};
% endfor

    // Copy all fluid-side gradients across to wall-side gradients
    ${pyfr.expand('bc_common_grad_copy', 'ub', 'nl', 'grad_ul', 'grad_ur')};

    // Remove the normal component of each species gradient
    fpdtype_t Ydotn;
% for n in range(ns):
    Ydotn = ${' + '.join(f'(grad_ul[{d}][{n}] - Ygi[{n}]*rho_{"xyz"[d]})'
                         f'*nl[{d}]' for d in range(ndims))};
% for d in range(ndims):
    grad_ur[${d}][${n}] -= Ydotn*nl[${d}];
% endfor
% endfor
</%pyfr:macro>
