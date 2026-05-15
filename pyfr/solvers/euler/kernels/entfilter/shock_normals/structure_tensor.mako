<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Structure-tensor shock normal from the low-order-projected polynomial
// gradient.  Uses the same per-upt physical gradient as volume_grad, but
// aggregates via the symmetric 2-tensor
//
//     S = sum_q  w_q * (grad phi_low)(x_q) (x) (grad phi_low)(x_q)
//
// and reports the dominant eigenvector of S.  Compared to the weighted-mean
// approach this is invariant under sign flips of the per-upt gradient
// (Gibbs lobes pointing in opposite directions both contribute positively
// to S rather than cancelling) and the dominant eigenvalue is itself the
// physical "shock-ness" weight.
//
// The low-order projection (shock-dir-p in cfg) is already baked into grad_op
// at the Python level, so the gradient values computed here are already the
// Gibbs-suppressed gradient.  The structure-tensor magnitude `n_mag` is
// sqrt(lambda_max), the dominant singular value of the gradient cloud.
//
// 2D only currently; 3D would require a 3x3 symmetric eigensolver.

<%pyfr:macro name='shock_normal_structure_tensor'
             params='u, grad_op, smats_upts, rcpdjac_upts, n_phys, n_mag'>
% if ndims != 2:
<% raise NotImplementedError('structure_tensor detector is 2D only currently') %>
% endif

    // Build per-upt scalar field phi (same field-choice convention as
    // volume_grad: shared vg_field cfg key).
    fpdtype_t phi[${nupts}];
    % for i in range(nupts):
    % if vg_field == 'density':
    phi[${i}] = u[${i}][0];
    % elif vg_field == 'pressure':
    {
        fpdtype_t rho = u[${i}][0];
        fpdtype_t E = u[${i}][${nvars - 1}];
        fpdtype_t ke = 0.5*(${' + '.join(
            f'u[{i}][{k}]*u[{i}][{k}]' for k in range(1, ndims + 1))})/rho;
        phi[${i}] = ${c['gamma'] - 1}*(E - ke);
    }
    % elif vg_field == 'entropy':
    {
        fpdtype_t rho = u[${i}][0];
        fpdtype_t E = u[${i}][${nvars - 1}];
        fpdtype_t ke = 0.5*(${' + '.join(
            f'u[{i}][{k}]*u[{i}][{k}]' for k in range(1, ndims + 1))})/rho;
        fpdtype_t p = ${c['gamma'] - 1}*(E - ke);
        phi[${i}] = (rho > 0 && p > 0)
                    ? p*pow(1.0/rho, ${c['gamma']}) : 0.0;
    }
    % else:
    <% raise ValueError(f"Unknown vg_field: {vg_field!r}") %>
    % endif
    % endfor

    // Reference-space gradient at every upt (grad_op is already
    // m4 @ projection, so this is the low-order-projected gradient)
    fpdtype_t grad_ref[${nupts}][${ndims}];
    % for q, d in pyfr.ndrange(nupts, ndims):
    grad_ref[${q}][${d}] = ${pyfr.dot(f'grad_op[{d*nupts + q}][{{k}}]',
                                       'phi[{k}]', k=nupts)};
    % endfor

    // Convert to physical gradient via J^{-T} = smats^T * rcpdjac
    fpdtype_t grad_phys[${nupts}][${ndims}];
    % for q in range(nupts):
    % for j in range(ndims):
    grad_phys[${q}][${j}] = (${' + '.join(
        f'smats_upts[{q}][{k*ndims + j}]*grad_ref[{q}][{k}]'
        for k in range(ndims))})*rcpdjac_upts[${q}];
    % endfor
    % endfor

    // Accumulate structure tensor S = sum_q w_q * g_q (x) g_q
    // and average gradient g_avg = sum_q w_q * g_q (for sign disambiguation)
    fpdtype_t S00 = 0.0, S01 = 0.0, S11 = 0.0;
    fpdtype_t g_avg0 = 0.0, g_avg1 = 0.0;
    % for q in range(nupts):
    {
        fpdtype_t w_q = (${upts_wts[q]})/rcpdjac_upts[${q}];
        fpdtype_t gx = grad_phys[${q}][0];
        fpdtype_t gy = grad_phys[${q}][1];
        S00 += w_q*gx*gx;
        S01 += w_q*gx*gy;
        S11 += w_q*gy*gy;
        g_avg0 += w_q*gx;
        g_avg1 += w_q*gy;
    }
    % endfor

    // 2x2 symmetric eigenvalue: S = [[a, b], [b, c]]
    //   lambda_max = (a+c)/2 + sqrt(((a-c)/2)^2 + b^2)
    fpdtype_t tr_half = 0.5*(S00 + S11);
    fpdtype_t diff_half = 0.5*(S00 - S11);
    fpdtype_t disc = sqrt(diff_half*diff_half + S01*S01);
    fpdtype_t lam_max = tr_half + disc;

    n_mag = sqrt(fmax(lam_max, 0.0));

    if (n_mag > ${shock_normal_eps})
    {
        // Eigenvector for lam_max.  Use (b, lam_max - a) or (lam_max - c, b)
        // whichever has the larger denominator for numerical stability.
        fpdtype_t ex, ey;
        if (fabs(S01) > ${shock_normal_eps})
        {
            ex = lam_max - S11;
            ey = S01;
        }
        else
        {
            // Diagonal: pick whichever entry is larger
            if (S00 > S11) { ex = 1.0; ey = 0.0; }
            else           { ex = 0.0; ey = 1.0; }
        }
        fpdtype_t inv_norm = 1.0/sqrt(ex*ex + ey*ey);
        ex *= inv_norm;
        ey *= inv_norm;

        // Sign: align with g_avg so n points toward high-phi region
        if (g_avg0*ex + g_avg1*ey < 0.0) { ex = -ex; ey = -ey; }

        n_phys[0] = ex;
        n_phys[1] = ey;
    }
    else
    {
        n_phys[0] = 0.0;
        n_phys[1] = 0.0;
    }
</%pyfr:macro>
