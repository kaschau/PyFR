<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.filters.iso_exponential'/>

// Pointwise admissibility limiter against the iso_exponential filter.
//
// Builds the modal coefficients, then walks every evaluation point shrinking
// the per-cell f via bisection whenever (rho, p, e) miss their bounds at
// that point.  After the scan f is the tightest value compatible with every
// constraint at every point; the filter is then applied in full at f and the
// mass defect from reference-space filtering is corrected against uavg.

<%pyfr:macro name='limiter_pointwise_bisect_lumped'
             params='u, uavg, entmin, f, vdm, invvdm, mean_wts'>
    // Compute modal basis
    fpdtype_t umodes[${nupts}][${nvars}];
    for (int uidx = 0; uidx < ${nupts}; uidx++)
    {
        for (int vidx = 0; vidx < ${nvars}; vidx++)
        {
            umodes[uidx][vidx] = ${pyfr.dot('invvdm[uidx][{k}]', 'u[{k}][vidx]', k=nupts)};
        }
    }

    // Setup filter (solve for f = exp(-zeta))
    fpdtype_t f_low, f_high, fnew;

    fpdtype_t d, p, e;

    // Compute f on a rolling basis per solution point
    fpdtype_t up[${order+1}][${nvars}];

    for (int uidx = 0; uidx < ${nefpts}; uidx++)
    {
        // Group nodal contributions by common filter factor
        % for pidx, vidx in pyfr.ndrange(order+1, nvars):
        up[${pidx}][${vidx}] = (${' + '.join(f'vdm[uidx][{k}]*umodes[{k}][{vidx}]'
                                               for k, dd in enumerate(ubdegs) if dd == pidx)});
        % endfor

        // Compute constraints with current minimum f value
        ${pyfr.expand('iso_exp_apply_single', 'up', 'f', 'd', 'p', 'e')};

        // Update f if constraints aren't satisfied
        if (d < ${d_min} || p < ${p_min} || e < entmin - ${e_tol})
        {
            // Set root-finding interval
            f_high = f;
            f_low = 0.0;

            // Iterate filter strength with bisection algorithm
            for (int iter = 0; iter < ${niters} && f_high - f_low > ${f_tol}; iter++)
            {
                // Compute new guess using bisection
                fnew = 0.5*(f_low + f_high);

                // Compute filtered state
                ${pyfr.expand('iso_exp_apply_single', 'up', 'fnew', 'd', 'p', 'e')};

                // Update brackets
                if (d < ${d_min} || p < ${p_min} || e < entmin - ${e_tol})
                    f_high = fnew;
                else
                    f_low = fnew;
            }

            // Set current minimum f as the bounds-preserving value
            f = f_low;
        }
    }

    // Filter full solution with bounds-preserving f value
    ${pyfr.expand('iso_exp_apply_full', 'umodes', 'vdm', 'u', 'f')};

    // Account for mass defect in reference-space filtering
    fpdtype_t duavg[${nvars}] = {0};
    % for vidx in range(nvars):
    duavg[${vidx}] = uavg[${vidx}] - ${pyfr.dot('mean_wts[{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
    % endfor

    %for uidx, vidx in pyfr.ndrange(nupts, nvars):
    u[${uidx}][${vidx}] += duavg[${vidx}];
    % endfor
</%pyfr:macro>
