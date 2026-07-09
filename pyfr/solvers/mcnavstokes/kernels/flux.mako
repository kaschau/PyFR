<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='viscous_flux_add'
             params='uin, grad_uin, rho, invrho, v, Y, T, cpmix, gammamix, mu, kappa, D, h, e_sp, fout'>
% for d in range(ndims):
    fpdtype_t rho_${'xyz'[d]} = ${' + '.join(f'grad_uin[{d}][{n}]'
                                             for n in range(ns))};
% endfor

    // Velocity derivatives
% for i, d in pyfr.ndrange(ndims, ndims):
    fpdtype_t v${i}_${'xyz'[d]} = invrho*(grad_uin[${d}][${ns + i}]
                                          - v[${i}]*rho_${'xyz'[d]});
% endfor

    // Temperature derivatives
    fpdtype_t E = uin[${nvars - 1}]*invrho;
    fpdtype_t rcpcv = gammamix/cpmix;
% for d in range(ndims):
    fpdtype_t eyy_${'xyz'[d]} = ${' + '.join(
        f'e_sp[{n}]*invrho*(grad_uin[{d}][{n}] - Y[{n}]*rho_{"xyz"[d]})'
        for n in range(ns))};
    fpdtype_t T_${'xyz'[d]} = rcpcv*(invrho*(grad_uin[${d}][${nvars - 1}]
        - E*rho_${'xyz'[d]})
        - (${' + '.join(f'v[{i}]*v{i}_{"xyz"[d]}' for i in range(ndims))})
        - eyy_${'xyz'[d]});
% endfor

    // Negated stress tensor
    fpdtype_t divv = ${' + '.join(f'v{i}_{"xyz"[i]}' for i in range(ndims))};
% for i, d in pyfr.ndrange(ndims, ndims):
% if i == d:
    fpdtype_t t_${i}${d} = -2*mu*(v${i}_${'xyz'[d]} - ${1.0/3.0}*divv);
% elif i < d:
    fpdtype_t t_${i}${d} = -mu*(v${i}_${'xyz'[d]} + v${d}_${'xyz'[i]});
% endif
% endfor

% for i, d in pyfr.ndrange(ndims, ndims):
    fout[${d}][${ns + i}] += t_${min(i, d)}${max(i, d)};
% endfor

    // Thermal diffusion
% for d in range(ndims):
    fout[${d}][${nvars - 1}] += ${' + '.join(
        f'v[{i}]*t_{min(i, d)}{max(i, d)}' for i in range(ndims))}
        + -kappa*T_${'xyz'[d]};
% endfor

    // Species diffusion
% for d in range(ndims):
    fpdtype_t Yd_${'xyz'[d]}[${ns}];
    fpdtype_t Vc_${'xyz'[d]} = 0.0;
% for n in range(ns):
    Yd_${'xyz'[d]}[${n}] = D[${n}]*(grad_uin[${d}][${n}]
                                    - Y[${n}]*rho_${'xyz'[d]});
    Vc_${'xyz'[d]} += Yd_${'xyz'[d]}[${n}];
% endfor
% endfor

    fpdtype_t J;
% for d, n in pyfr.ndrange(ndims, ns):
    J = -Yd_${'xyz'[d]}[${n}] + Y[${n}]*Vc_${'xyz'[d]};
    fout[${d}][${n}] += J;
    fout[${d}][${nvars - 1}] += h[${n}]*J;
% endfor
</%pyfr:macro>
