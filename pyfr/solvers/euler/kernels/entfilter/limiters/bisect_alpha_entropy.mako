<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Entropy limiter for the convex blend
//     u(α) = α·u_orig + (1-α)·u_target.
//
// Entropy s(α) = p(α)·ρ(α)^(-γ) is nonlinear in α, so we bisect.  For each
// admissibility-violating eval point we find the largest α such that
//     s(α) ≥ entmin − e_tol
// then take the global minimum across all eval points.
//
// Per-(α, eval-point) state evaluation is cheap because the blend is a
// linear combination of two pre-computed nodal arrays — no modal
// reconstruction needed.

<%pyfr:macro name='solve_alpha_entropy'
             params='u_orig_ext, u_target_ext, entmin, alpha_e'>
    alpha_e = 1.0;
    fpdtype_t gm1 = ${c['gamma'] - 1};
    fpdtype_t s_threshold = entmin - ${e_tol};

    for (int eidx = 0; eidx < ${nefpts}; eidx++)
    {
        // Compute s at α=1 (orig).  If admissible, no constraint from this pt.
        fpdtype_t rho_o = u_orig_ext[eidx][0];
        fpdtype_t E_o = u_orig_ext[eidx][${nvars - 1}];
        fpdtype_t m2_o = ${' + '.join(
            f'u_orig_ext[eidx][{k}]*u_orig_ext[eidx][{k}]'
            for k in range(1, ndims + 1))};
        fpdtype_t p_o = gm1*(E_o - 0.5*m2_o/rho_o);
        fpdtype_t s_o = (rho_o > 0 && p_o > 0)
                        ? p_o*pow(1.0/rho_o, ${c['gamma']})
                        : -${fpdtype_max};

        if (s_o >= s_threshold) continue;     // already admissible at α=1

        // Bisect on α in [0, 1].  s(α=0) = s_target; assume target admissible.
        fpdtype_t a_lo = 0.0;
        fpdtype_t a_hi = 1.0;
        fpdtype_t a_mid;
        for (int iter = 0; iter < ${niters} && a_hi - a_lo > ${f_tol}; iter++)
        {
            a_mid = 0.5*(a_lo + a_hi);

            // Evaluate state at α=a_mid by linear blend
            fpdtype_t rho_m = a_mid*rho_o + (1.0 - a_mid)*u_target_ext[eidx][0];
            fpdtype_t E_m = a_mid*E_o + (1.0 - a_mid)*u_target_ext[eidx][${nvars - 1}];
            fpdtype_t m2_m = ${' + '.join(
                f'(a_mid*u_orig_ext[eidx][{k}] + (1.0 - a_mid)*u_target_ext[eidx][{k}])*'
                f'(a_mid*u_orig_ext[eidx][{k}] + (1.0 - a_mid)*u_target_ext[eidx][{k}])'
                for k in range(1, ndims + 1))};
            fpdtype_t p_m = gm1*(E_m - 0.5*m2_m/rho_m);
            fpdtype_t s_m = (rho_m > 0 && p_m > 0)
                            ? p_m*pow(1.0/rho_m, ${c['gamma']})
                            : -${fpdtype_max};

            if (s_m < s_threshold) { a_hi = a_mid; }   // tighten (lower α)
            else                   { a_lo = a_mid; }   // relax (higher α)
        }

        alpha_e = fmin(alpha_e, a_lo);
    }
    alpha_e = fmax(0.0, fmin(1.0, alpha_e));
</%pyfr:macro>
