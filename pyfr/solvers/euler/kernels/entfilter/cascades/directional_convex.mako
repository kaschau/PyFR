<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.shared.minima'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.limiters.linearised_dpe'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.shock_normals.${shock_normal_detector}'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.filters.directional_convex'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.limiters.closed_form_density'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.limiters.closed_form_pressure'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.limiters.bisect_alpha_entropy'/>

// Directional convex-blend cascade.  When a cell fails admissibility:
//   1. Detect shock normal n_ref (reference space) via configured detector.
//   2. Build u_target via the directional convex filter at full strength
//      (α=0 endpoint) — the modes perpendicular to n_ref in modal-index
//      space are preserved; others are exponentially damped.
//   3. Interpolate u to all eval points (upts + fpts) to get u_orig_ext.
//   4. Solve sequentially:
//        - α_d via closed form on density (linear in α).  Apply blend.
//        - α_p via closed form on pressure (quadratic in α).  Apply blend.
//        - α_e via bisection on entropy.  Apply blend.
//   5. Write the final blended state back to u (upts portion).
//   6. Fall through to linearised dpe if any constraint still violated.

<%pyfr:kernel name='entropyfilter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='inout fpdtype_t[${str(nfaces)}]'
              ef_filter='out fpdtype_t[1]'
              vdm='in broadcast fpdtype_t[${str(nefpts)}][${str(nupts)}]'
              invvdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'
              mean_wts='in fpdtype_t[${str(nupts)}]'
              smats_upts='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              rcpdjac_upts='in fpdtype_t[${str(nupts)}]'
              grad_op='in broadcast fpdtype_t[${str(ndims*nupts)}][${str(nupts)}]'
              shock_normal='out fpdtype_t[${str(ndims)}]'
              shock_normal_mag='out fpdtype_t[1]'>
    fpdtype_t dmin, pmin, emin;
    fpdtype_t f = 1.0;

    // Compute minimum entropy from current and adjacent elements
    fpdtype_t entmin = ${fpdtype_max};
    for (int fidx = 0; fidx < ${nfaces}; fidx++) entmin = fmin(entmin, entmin_int[fidx]);

    // Check if solution is within bounds
    ${pyfr.expand('get_minima', 'u', 'm0', 'dmin', 'pmin', 'emin')};

    // Compute mean quantities (used by linearised fallback)
    fpdtype_t uavg[${nvars}];
    % for vidx in range(nvars):
    uavg[${vidx}] = ${pyfr.dot('mean_wts[{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
    % endfor

    // Default outputs: zero for smooth cells
    fpdtype_t n_phys[${ndims}];
    fpdtype_t n_mag = 0.0;
    % for d in range(ndims):
    n_phys[${d}] = 0.0;
    % endfor

    if (dmin < ${d_min} || pmin < ${p_min} || emin < entmin - ${e_tol})
    {
        // ----------------------------------------------------------------
        // 1. Shock normal detection (ref-space n_ref + physical n_phys)
        // ----------------------------------------------------------------
        fpdtype_t n_ref[${ndims}];
        % for d in range(ndims):
        n_ref[${d}] = 0.0;
        % endfor

        % if shock_normal_detector == 'face_grad':
        ${pyfr.expand('shock_normal_face_grad', 'u', 'm0', 'n_ref', 'n_mag')};
        % elif shock_normal_detector == 'volume_grad_density':
        // volume_grad outputs physical-space directly; convert back to ref
        // by applying J (∝ smats^-T) — or equivalently, recompute in ref
        // space.  For now just emit a warning; user should use face_grad.
        <% raise NotImplementedError("Use shock-normal = face_grad with directional_convex cascade for now") %>
        % elif shock_normal_detector == 'face_entmin':
        ${pyfr.expand('shock_normal_face_entmin', 'entmin_int', 'n_ref', 'n_mag')};
        % else:
        <% raise ValueError(f"Unknown shock-normal detector: {shock_normal_detector!r}") %>
        % endif

        // Convert n_ref → physical for export visualisation (smats^T*n_ref)
        if (n_mag > ${shock_normal_eps})
        {
            fpdtype_t n_phys_raw[${ndims}];
            % for k in range(ndims):
            n_phys_raw[${k}] = ${' + '.join(
                f'smats_upts[0][{i*ndims + k}]*n_ref[{i}]'
                for i in range(ndims))};
            % endfor
            fpdtype_t inv_mag = 1.0/sqrt(${' + '.join(
                f'n_phys_raw[{k}]*n_phys_raw[{k}]' for k in range(ndims))});
            % for k in range(ndims):
            n_phys[${k}] = n_phys_raw[${k}]*inv_mag;
            % endfor
        }

        // ----------------------------------------------------------------
        // 2. Build u_target via directional convex filter at all eval points
        // ----------------------------------------------------------------
        fpdtype_t u_target_ext[${nefpts}][${nvars}];
        if (n_mag > ${shock_normal_eps})
        {
            ${pyfr.expand('build_u_target_directional_convex',
                           'u', 'invvdm', 'vdm', 'n_ref', 'u_target_ext')};
        }
        else
        {
            // Detector found nothing usable; fall back to cell-mean target
            // (this is the iso_linear_to_mean target — degenerate of the
            // directional filter with no preferred direction).
            for (int eidx = 0; eidx < ${nefpts}; eidx++)
            {
                % for vidx in range(nvars):
                u_target_ext[eidx][${vidx}] = uavg[${vidx}];
                % endfor
            }
        }

        // ----------------------------------------------------------------
        // 3. Build u_orig_ext: u at upts + interpolated u at fpts via m0
        // ----------------------------------------------------------------
        fpdtype_t u_orig_ext[${nefpts}][${nvars}];
        for (int uidx = 0; uidx < ${nupts}; uidx++)
        {
            % for vidx in range(nvars):
            u_orig_ext[uidx][${vidx}] = u[uidx][${vidx}];
            % endfor
        }
        % if not fpts_in_upts:
        for (int fidx = 0; fidx < ${nfpts}; fidx++)
        {
            % for vidx in range(nvars):
            u_orig_ext[${nupts} + fidx][${vidx}] = ${pyfr.dot('m0[fidx][{k}]',
                                                              f'u[{{k}}][{vidx}]', k=nupts)};
            % endfor
        }
        % endif

        // ----------------------------------------------------------------
        // 4. Sequential cascade: α_d, α_p, α_e
        // ----------------------------------------------------------------
        // Density (closed form)
        fpdtype_t alpha_d;
        ${pyfr.expand('solve_alpha_density', 'u_orig_ext', 'u_target_ext', 'alpha_d')};

        // Apply blend at α_d to u_orig_ext in-place
        if (alpha_d < 1.0)
        {
            for (int eidx = 0; eidx < ${nefpts}; eidx++)
            {
                % for vidx in range(nvars):
                u_orig_ext[eidx][${vidx}] = alpha_d*u_orig_ext[eidx][${vidx}]
                    + (1.0 - alpha_d)*u_target_ext[eidx][${vidx}];
                % endfor
            }
        }

        // Pressure (closed form)
        fpdtype_t alpha_p;
        ${pyfr.expand('solve_alpha_pressure', 'u_orig_ext', 'u_target_ext', 'alpha_p')};

        if (alpha_p < 1.0)
        {
            for (int eidx = 0; eidx < ${nefpts}; eidx++)
            {
                % for vidx in range(nvars):
                u_orig_ext[eidx][${vidx}] = alpha_p*u_orig_ext[eidx][${vidx}]
                    + (1.0 - alpha_p)*u_target_ext[eidx][${vidx}];
                % endfor
            }
        }

        // Entropy (bisection)
        fpdtype_t alpha_e;
        ${pyfr.expand('solve_alpha_entropy',
                       'u_orig_ext', 'u_target_ext', 'entmin', 'alpha_e')};

        if (alpha_e < 1.0)
        {
            for (int eidx = 0; eidx < ${nefpts}; eidx++)
            {
                % for vidx in range(nvars):
                u_orig_ext[eidx][${vidx}] = alpha_e*u_orig_ext[eidx][${vidx}]
                    + (1.0 - alpha_e)*u_target_ext[eidx][${vidx}];
                % endfor
            }
        }

        // ----------------------------------------------------------------
        // 5. Write final state back to u (upts portion)
        // ----------------------------------------------------------------
        for (int uidx = 0; uidx < ${nupts}; uidx++)
        {
            % for vidx in range(nvars):
            u[uidx][${vidx}] = u_orig_ext[uidx][${vidx}];
            % endfor
        }

        // Track combined filter strength (1 = no filter, 0 = max).
        // Approximate as the product of (1-effective-α) -- log-additive blending.
        f = fmin(f, alpha_d*alpha_p*alpha_e);

        // Refresh minima for fallback check
        ${pyfr.expand('get_minima', 'u', 'm0', 'dmin', 'pmin', 'emin')};
    }

    // Fallback: linearised dpe limiter if anything still out of bounds
    if (dmin < ${d_min} || pmin < ${p_min} || emin < entmin - ${e_tol})
    {
        ${pyfr.expand('limiter_linearised_dpe',
                       'u', 'uavg', 'entmin', 'dmin', 'pmin', 'emin', 'f', 'm0')};
    }

    // Set new minimum entropy within element for next stage
    % for fidx in range(nfaces):
    entmin_int[${fidx}] = emin;
    % endfor

    // Output the filter strength
    ef_filter[0] = f;

    // Output the shock normal (physical-space, unit-length where active)
    % for d in range(ndims):
    shock_normal[${d}] = n_phys[${d}];
    % endfor
    shock_normal_mag[0] = n_mag;
</%pyfr:kernel>
