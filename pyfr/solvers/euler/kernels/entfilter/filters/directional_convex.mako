<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Directional convex-blend filter: builds the fully-damped target state
// u_target = vdm @ (d ⊙ invvdm @ u), where
//     d_k = exp(-ζ_max · Q_k)
//     Q_k = (i_k · n_ξ + j_k · n_η)²
//
// Modes perpendicular to n_ref in modal-index space have Q_k = 0 and are
// preserved exactly; modes with Q_k > 0 are exponentially damped.  The
// constant mode (i=j=0) is always preserved.
//
// ζ_max is chosen per-cell so the worst-damped mode is reduced by 1/ε:
//     ζ_max = ln(1/ε) / Q_max
// with Q_max = max_k Q_k and ε from cfg (default 1e-3).
//
// The blend u(α) = α·u_orig + (1-α)·u_target is linear in α — caller can
// either compute u_target once and blend, or evaluate u(α) directly at
// each eval point.

<%pyfr:macro name='build_u_target_directional_convex'
             params='u, invvdm, vdm, n_ref, u_target'>
    // 1. Compute Q_k = (i·n_x + j·n_y)² for each mode (template-time
    //    expansion using mode_ij tplarg)
    fpdtype_t Q[${nupts}];
    fpdtype_t Q_max = 0.0;
    % for k, (i, j) in enumerate(mode_ij):
    % if ndims == 2:
    {
        fpdtype_t dot_k = ${i}*n_ref[0] + ${j}*n_ref[1];
        Q[${k}] = dot_k*dot_k;
        if (Q[${k}] > Q_max) Q_max = Q[${k}];
    }
    % elif ndims == 3:
    {
        fpdtype_t dot_k = ${i}*n_ref[0] + ${j}*n_ref[1] + ${mode_ij[k][2]}*n_ref[2];
        Q[${k}] = dot_k*dot_k;
        if (Q[${k}] > Q_max) Q_max = Q[${k}];
    }
    % endif
    % endfor

    // 2. Compute ζ_max from Q_max and chosen residual ε
    fpdtype_t zeta_max = log(1.0/${dir_eps})/fmax(Q_max, 1e-12);

    // 3. Compute per-mode damping d_k = exp(-ζ_max · Q_k)
    fpdtype_t d_mode[${nupts}];
    % for k in range(nupts):
    d_mode[${k}] = exp(-zeta_max*Q[${k}]);
    % endfor

    // 4. Modal coefficients: û = invvdm @ u
    fpdtype_t umodes[${nupts}][${nvars}];
    for (int uidx = 0; uidx < ${nupts}; uidx++)
    {
        for (int vidx = 0; vidx < ${nvars}; vidx++)
        {
            umodes[uidx][vidx] = ${pyfr.dot('invvdm[uidx][{k}]',
                                            'u[{k}][vidx]', k=nupts)};
        }
    }

    // 5. Damp each mode and recover nodal at all eval points
    //    (upts + fpts; nefpts total).  vdm has shape (nefpts, nupts).
    for (int eidx = 0; eidx < ${nefpts}; eidx++)
    {
        for (int vidx = 0; vidx < ${nvars}; vidx++)
        {
            fpdtype_t tmp = 0.0;
            for (int k = 0; k < ${nupts}; k++)
            {
                tmp += vdm[eidx][k]*d_mode[k]*umodes[k][vidx];
            }
            u_target[eidx][vidx] = tmp;
        }
    }
</%pyfr:macro>
