<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.shared.minima'/>
<%include file='pyfr.solvers.euler.kernels.entfilter.filters.iso_linear_to_mean'/>

// Sequential rho -> p -> e limiter against the iso_linear_to_mean filter.
//
// For each bound that's violated, compute the closed-form blend coefficient
// alpha toward the cell mean, apply, and recompute minima.  Density only
// touches u[*][0] (the other variables are insensitive to the rho bound);
// pressure and entropy touch all variables.

<%pyfr:macro name='limiter_linearised_dpe'
             params='u, uavg, entmin, dmin, pmin, emin, f, m0'>
    fpdtype_t davg, pavg, eavg;
    ${pyfr.expand('compute_entropy', 'uavg', 'davg', 'pavg', 'eavg')};

    // Apply density, pressure, and entropy limiting sequentially
    fpdtype_t alpha;
    % for (fvar, bound, blend) in [('d', d_min, 'iso_lin_blend_density'), ('p', p_min, 'iso_lin_blend_full'), ('e', f'entmin - {e_tol}', 'iso_lin_blend_full')]:
    if (${fvar}min < ${bound})
    {
        alpha = (${fvar}min - (${bound}))/(${fvar}min - ${fvar}avg);
        alpha = fmin(fmax(alpha, 0.0), 1.0);
        f = fmin(f, 1.0 - alpha);

        ${pyfr.expand(blend, 'u', 'uavg', 'alpha')};

        ${pyfr.expand('get_minima', 'u', 'm0', 'dmin', 'pmin', 'emin')};
    }
    % endfor
</%pyfr:macro>
