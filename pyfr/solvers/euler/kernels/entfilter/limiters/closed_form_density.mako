<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Closed-form density limiter for the convex blend
//     u(α) = α·u_orig + (1-α)·u_target.
//
// ρ(α) is linear in α at each eval point.  Constraint ρ(α) ≥ d_min:
//
//     α·(ρ_orig − ρ_target) ≥ d_min − ρ_target
//
// For the typical filter case (ρ_orig < d_min ≤ ρ_target), this gives
// the most-restrictive α at each violating eval point as
//
//     α_q = (d_min − ρ_target_q) / (ρ_orig_q − ρ_target_q)
//
// (denominator and numerator both negative → α_q ∈ (0, 1)).  The cell-
// wide limiter takes the minimum across violating points.  Caller is
// responsible for applying the resulting blend.

<%pyfr:macro name='solve_alpha_density'
             params='u_orig_ext, u_target_ext, alpha_d'>
    alpha_d = 1.0;
    for (int eidx = 0; eidx < ${nefpts}; eidx++)
    {
        fpdtype_t rho_o = u_orig_ext[eidx][0];
        if (rho_o < ${d_min})
        {
            fpdtype_t rho_t = u_target_ext[eidx][0];
            // α_q = (d_min − ρ_t)/(ρ_o − ρ_t); valid when ρ_o ≠ ρ_t.
            // For our setup ρ_o < d_min ≤ ρ_t so denominator < 0.
            fpdtype_t denom = rho_o - rho_t;
            if (fabs(denom) > 1e-30)
            {
                fpdtype_t alpha_q = (${d_min} - rho_t)/denom;
                alpha_d = fmin(alpha_d, alpha_q);
            }
            else
            {
                // Degenerate: ρ_t = ρ_o < d_min, target also inadmissible
                alpha_d = 0.0;
            }
        }
    }
    alpha_d = fmax(0.0, fmin(1.0, alpha_d));
</%pyfr:macro>
