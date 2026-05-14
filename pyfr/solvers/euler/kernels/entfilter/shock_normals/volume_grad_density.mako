<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Volume-weighted shock normal from the polynomial gradient of a chosen
// scalar field (density or entropy, selected via the `vg_field` tplarg).
//
// At each upt q, compute the physical gradient of phi:
//     grad_phys[q] = J^{-T} * grad_ref[q]
//                  = smats^T * rcpdjac * grad_ref[q]
// where grad_ref[q][d] = sum_i grad_op[d*nupts + q][i] * phi[i].
//
// Then form the steepness-weighted volume average:
//     grad_avg[j] = sum_q (w_q * s_q^p * grad_phys[q][j]) / sum_q (w_q * s_q^p)
// with w_q = upts_wts[q] * det(J)[q] = upts_wts[q] / rcpdjac_upts[q] and
// s_q = |grad_phys[q]|.  weight_power p emphasises steep upts; for p=0 it
// collapses to a flat volume-weighted mean gradient, p=1 weights by gradient
// magnitude, p=2 emphasises the steepest upts more aggressively.
//
// Output n_phys is the unit physical-space direction of grad_avg; n_mag is the
// magnitude of the weighted-mean gradient (useful as a shock-strength proxy).

<%pyfr:macro name='shock_normal_volume_grad_density'
             params='u, grad_op, smats_upts, rcpdjac_upts, n_phys, n_mag'>
    // Build per-upt scalar field phi (density or entropy)
    fpdtype_t phi[${nupts}];
    % for i in range(nupts):
    % if vg_field == 'density':
    phi[${i}] = u[${i}][0];
    % elif vg_field == 'entropy':
    {
        fpdtype_t rho = u[${i}][0];
        fpdtype_t rcprho = 1.0/rho;
        fpdtype_t E = u[${i}][${nvars - 1}];
        fpdtype_t p = ${c['gamma'] - 1}*(E - 0.5*rcprho*(${' + '.join(
            f'u[{i}][{k}]*u[{i}][{k}]' for k in range(1, ndims + 1))}));
        phi[${i}] = (rho > 0 && p > 0) ? p*pow(rcprho, ${c['gamma']}) : 0.0;
    }
    % else:
    <% raise ValueError(f"Unknown vg_field: {vg_field!r}") %>
    % endif
    % endfor

    // 1. Reference-space gradient of phi at every upt
    fpdtype_t grad_ref[${nupts}][${ndims}];
    % for q, d in pyfr.ndrange(nupts, ndims):
    grad_ref[${q}][${d}] = ${pyfr.dot(f'grad_op[{d*nupts + q}][{{k}}]',
                                       'phi[{k}]', k=nupts)};
    % endfor

    // 2. Convert to physical gradient at every upt via J^{-T} = smats^T * rcpdjac
    fpdtype_t grad_phys[${nupts}][${ndims}];
    fpdtype_t s[${nupts}];
    % for q in range(nupts):
    % for j in range(ndims):
    grad_phys[${q}][${j}] = (${' + '.join(
        f'smats_upts[{q}][{k*ndims + j}]*grad_ref[{q}][{k}]'
        for k in range(ndims))})*rcpdjac_upts[${q}];
    % endfor
    s[${q}] = sqrt(${' + '.join(
        f'grad_phys[{q}][{j}]*grad_phys[{q}][{j}]' for j in range(ndims))});
    % endfor

    // 3. Steepness-weighted volume average
    fpdtype_t total_w = 0.0;
    fpdtype_t grad_avg[${ndims}];
    % for j in range(ndims):
    grad_avg[${j}] = 0.0;
    % endfor

    % for q in range(nupts):
    {
        fpdtype_t w_q = (${upts_wts[q]})/rcpdjac_upts[${q}]
                        % if vg_weight_power == 0:
                        ;
                        % elif vg_weight_power == 1:
                        *s[${q}];
                        % else:
                        *pow(s[${q}], ${float(vg_weight_power)});
                        % endif
        total_w += w_q;
        % for j in range(ndims):
        grad_avg[${j}] += w_q*grad_phys[${q}][${j}];
        % endfor
    }
    % endfor

    // 4. Normalise
    if (total_w > ${shock_normal_eps})
    {
        fpdtype_t inv_w = 1.0/total_w;
        % for j in range(ndims):
        grad_avg[${j}] *= inv_w;
        % endfor

        n_mag = sqrt(${' + '.join(
            f'grad_avg[{j}]*grad_avg[{j}]' for j in range(ndims))});

        if (n_mag > ${shock_normal_eps})
        {
            fpdtype_t inv_mag = 1.0/n_mag;
            % for j in range(ndims):
            n_phys[${j}] = grad_avg[${j}]*inv_mag;
            % endfor
        }
        else
        {
            % for j in range(ndims):
            n_phys[${j}] = 0.0;
            % endfor
        }
    }
    else
    {
        n_mag = 0.0;
        % for j in range(ndims):
        n_phys[${j}] = 0.0;
        % endfor
    }
</%pyfr:macro>
