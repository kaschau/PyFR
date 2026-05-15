<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Closed-form pressure limiter for the convex blend
//     u(α) = α·u_orig + (1-α)·u_target.
//
// Pressure p = (γ-1)·(E - 0.5·|m|²/ρ).  The constraint p(α) ≥ p_min,
// multiplied through by ρ(α) > 0, becomes a quadratic in α:
//
//     F(α) = A + B·α + C·α² ≥ 0
//
// where (with δX = X_orig − X_target for each conserved variable X):
//     A = ρ_t·(p_t − p_min)
//     B = (γ-1)·[E_t·δρ + δE·ρ_t − Σ_i m_t_i·δm_i] − p_min·δρ
//     C = (γ-1)·δE·δρ − 0.5·(γ-1)·Σ_i δm_i²
//
// For the typical filter case (p_orig < p_min, p_target ≥ p_min):
//   F(0) = ρ_t·(p_t − p_min) ≥ 0
//   F(1) = ρ_o·(p_orig − p_min) < 0
// so F crosses zero somewhere in (0, 1].  Solve via quadratic formula
// and pick the smallest non-negative root in [0, 1].

<%pyfr:macro name='solve_alpha_pressure'
             params='u_orig_ext, u_target_ext, alpha_p'>
    alpha_p = 1.0;
    fpdtype_t gm1 = ${c['gamma'] - 1};
    for (int eidx = 0; eidx < ${nefpts}; eidx++)
    {
        fpdtype_t rho_o = u_orig_ext[eidx][0];
        fpdtype_t rho_t = u_target_ext[eidx][0];
        fpdtype_t E_o = u_orig_ext[eidx][${nvars - 1}];
        fpdtype_t E_t = u_target_ext[eidx][${nvars - 1}];

        // |m|² for orig and target
        fpdtype_t m2_o = ${' + '.join(
            f'u_orig_ext[eidx][{k}]*u_orig_ext[eidx][{k}]'
            for k in range(1, ndims + 1))};
        fpdtype_t m2_t = ${' + '.join(
            f'u_target_ext[eidx][{k}]*u_target_ext[eidx][{k}]'
            for k in range(1, ndims + 1))};

        // p at α=1 (orig) and α=0 (target)
        fpdtype_t p_o = gm1*(E_o - 0.5*m2_o/rho_o);
        if (p_o >= ${p_min}) continue;          // already admissible at α=1

        fpdtype_t p_t = gm1*(E_t - 0.5*m2_t/rho_t);

        // Deltas
        fpdtype_t drho = rho_o - rho_t;
        fpdtype_t dE = E_o - E_t;
        // Σ_i m_t_i·δm_i  and  Σ_i δm_i²
        fpdtype_t dot_mt_dm = ${' + '.join(
            f'u_target_ext[eidx][{k}]*(u_orig_ext[eidx][{k}]-u_target_ext[eidx][{k}])'
            for k in range(1, ndims + 1))};
        fpdtype_t sum_dm2 = ${' + '.join(
            f'(u_orig_ext[eidx][{k}]-u_target_ext[eidx][{k}])'
            f'*(u_orig_ext[eidx][{k}]-u_target_ext[eidx][{k}])'
            for k in range(1, ndims + 1))};

        fpdtype_t A = rho_t*(p_t - ${p_min});
        fpdtype_t B = gm1*(E_t*drho + dE*rho_t - dot_mt_dm) - ${p_min}*drho;
        fpdtype_t C = gm1*(dE*drho - 0.5*sum_dm2);

        // Solve C·α² + B·α + A = 0 for smallest non-negative root in [0, 1].
        fpdtype_t alpha_q;
        if (fabs(C) < 1e-30)
        {
            // Linear case: B·α + A = 0  →  α = -A/B
            if (fabs(B) < 1e-30) { alpha_q = 0.0; }
            else { alpha_q = -A/B; }
        }
        else
        {
            fpdtype_t disc = B*B - 4.0*A*C;
            if (disc < 0.0) { alpha_q = 0.0; }     // shouldn't happen
            else
            {
                fpdtype_t sd = sqrt(disc);
                fpdtype_t r1 = (-B + sd)/(2.0*C);
                fpdtype_t r2 = (-B - sd)/(2.0*C);
                // Want the smallest root in [0, 1].  Since F(0)≥0, F(1)<0
                // the relevant root is whichever lies in (0, 1].
                fpdtype_t lo = fmin(r1, r2);
                fpdtype_t hi = fmax(r1, r2);
                alpha_q = (lo >= 0.0 && lo <= 1.0) ? lo
                         : ((hi >= 0.0 && hi <= 1.0) ? hi : 0.0);
            }
        }
        alpha_q = fmax(0.0, fmin(1.0, alpha_q));
        alpha_p = fmin(alpha_p, alpha_q);
    }
</%pyfr:macro>
