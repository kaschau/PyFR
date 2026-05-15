<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.shared.minima'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.limiters.pointwise_bisect_lumped'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.limiters.linearised_dpe'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.shock_normals.${shock_normal_detector}'/>

// Today's "nonlinear" cascade:  pointwise bisection on the iso_exponential
// filter, with the linearised dpe limiter as a final fallback if anything is
// still out of bounds after the bisection step.
//
// Side-effect: when admissibility fails, compute the shock normal using the
// detector named by `shock_normal_detector` (cfg key shock-normal) and write
// to shock_normal/shock_normal_mag (physical-space unit vector + magnitude).
// Smooth cells get zeros.

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

    // Compute mean quantities using per-element weights
    fpdtype_t uavg[${nvars}];
    % for vidx in range(nvars):
    uavg[${vidx}] = ${pyfr.dot('mean_wts[{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
    % endfor

    // Default shock-normal output: zero for smooth cells.
    fpdtype_t n_phys[${ndims}];
    fpdtype_t n_mag = 0.0;
    % for d in range(ndims):
    n_phys[${d}] = 0.0;
    % endfor

    // Bisection limiter against the iso_exponential filter
    if (dmin < ${d_min} || pmin < ${p_min} || emin < entmin - ${e_tol})
    {
        // Compute shock normal as a side-effect of admissibility failure.
        // Each detector exposes a macro `compute_shock_normal` that writes
        // the physical-space unit normal into n_phys and the magnitude proxy
        // into n_mag, with whatever inputs its strategy needs.
        % if shock_normal_detector == 'face_entmin':
        {
            // face_entmin builds a reference-space normal from per-face
            // entmin via the divergence theorem; convert to physical via
            // smats^T (using upt 0 -- constant per cell for linear elements).
            fpdtype_t n_ref[${ndims}];
            ${pyfr.expand('shock_normal_face_entmin',
                           'entmin_int', 'n_ref', 'n_mag')};
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
        }
        % elif shock_normal_detector == 'volume_grad_density':
        ${pyfr.expand('shock_normal_volume_grad_density',
                       'u', 'grad_op', 'smats_upts', 'rcpdjac_upts',
                       'n_phys', 'n_mag')};
        % elif shock_normal_detector == 'face_grad':
        {
            // face_grad outputs reference-space; convert to physical via
            // smats^T (using upt 0 -- constant per cell for linear elements).
            fpdtype_t n_ref[${ndims}];
            ${pyfr.expand('shock_normal_face_grad', 'u', 'm0',
                           'n_ref', 'n_mag')};
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
        }
        % elif shock_normal_detector == 'structure_tensor':
        ${pyfr.expand('shock_normal_structure_tensor',
                       'u', 'grad_op', 'smats_upts', 'rcpdjac_upts',
                       'n_phys', 'n_mag')};
        % else:
        <% raise ValueError(f"Unknown shock-normal detector: {shock_normal_detector!r}") %>
        % endif

        ${pyfr.expand('limiter_pointwise_bisect_lumped',
                       'u', 'uavg', 'entmin', 'f', 'vdm', 'invvdm', 'mean_wts')};

        // Calculate minimum entropy from filtered solution
        ${pyfr.expand('get_minima', 'u', 'm0', 'dmin', 'pmin', 'emin')};
    }

    // Linearised limiter against the iso_linear_to_mean filter as fallback
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

    // Output the shock normal (physical-space, unit-length where active;
    // zero in smooth cells) and its magnitude (proxy for shock strength).
    % for d in range(ndims):
    shock_normal[${d}] = n_phys[${d}];
    % endfor
    shock_normal_mag[0] = n_mag;
</%pyfr:kernel>
