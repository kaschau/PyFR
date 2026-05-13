<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.shared.minima'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.limiters.pointwise_bisect_lumped'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.limiters.linearised_dpe'/>

// Today's "nonlinear" cascade:  pointwise bisection on the iso_exponential
// filter, with the linearised dpe limiter as a final fallback if anything is
// still out of bounds after the bisection step.

<%pyfr:kernel name='entropyfilter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='inout fpdtype_t[${str(nfaces)}]'
              ef_filter='out fpdtype_t[1]'
              vdm='in broadcast fpdtype_t[${str(nefpts)}][${str(nupts)}]'
              invvdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'
              mean_wts='in fpdtype_t[${str(nupts)}]'>
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

    // Bisection limiter against the iso_exponential filter
    if (dmin < ${d_min} || pmin < ${p_min} || emin < entmin - ${e_tol})
    {
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
</%pyfr:kernel>
