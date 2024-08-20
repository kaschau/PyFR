<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.entropy'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='get_min_rhoY' params='u, q, rhoYmin'>

    rhoYmin = ${fpdtype_max};
    % for n in range(ns - 1):
    rhoYmin = fmin(rhoYmin, u[${Yix + n}]);
    % endfor
    // Check ns species
    rhoYmin = fmin(rhoYmin, q[${Yix + ns - 1}]*u[0]);

</%pyfr:macro>

<%pyfr:macro name='get_minima' params='u, m0, rhomin, rhoYmin, pmin, emin'>

    rhomin = ${fpdtype_max};
    rhoYmin = ${fpdtype_max};
    pmin = ${fpdtype_max};
    emin = ${fpdtype_max};

    for (int i = 0; i < ${nupts}; i++)
    {
    fpdtype_t ui[${nvars}];
    % for j in range(nvars):
        ui[${j}] = u[i][${j}];
    % endfor

        // Compute thermodynamic properties
        fpdtype_t qi[${nvars+1}];
        fpdtype_t qhi[${3+ns}];
        ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

        fpdtype_t Ymintemp;
        ${pyfr.expand('get_min_Y', 'qi', 'Ymintemp')};
        fpdtype_t e;
        ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};
        fpdtype_t rhoYmintemp;
        ${pyfr.expand('get_min_rhoY', 'ui', 'qi', 'rhoYmintemp')};

        rhomin = fmin(rhomin, ui[0]);
        rhoYmin = fmin(rhoYmin, rhoYmintemp);
        pmin = fmin(pmin, qi[0]);
        emin = fmin(emin, e);
    }

    % if not fpts_in_upts:
    fpdtype_t uf[${nvars}];
    for (int fidx = 0; fidx < ${nfpts}; fidx++)
    {
        % for vidx in range(nvars):
        uf[${vidx}] = ${pyfr.dot('m0[fidx][{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
        % endfor

        // Compute thermodynamic properties
        fpdtype_t qf[${nvars+1}];
        fpdtype_t qhf[${3+ns}];
        ${pyfr.expand('stateFrom-cons', 'uf', 'qf', 'qhf')};

        fpdtype_t e;
        ${pyfr.expand('compute_entropy', 'uf', 'qf', 'e')};
        fpdtype_t rhoYmintemp;
        ${pyfr.expand('get_min_rhoY', 'uf', 'qf', 'rhoYmintemp')};

        rhomin = fmin(rhomin, ui[0]);
        rhoYmin = fmin(rhoYmin, rhoYmintemp);
        pmin = fmin(pmin, qi[0]);
        emin = fmin(emin, e);
    }
    % endif
</%pyfr:macro>

<%pyfr:macro name='apply_filter_full' params='umodes, vdm, uf, f'>
    // Precompute filter factors per basis degree
    fpdtype_t ffac[${order + 1}];
    fpdtype_t v = ffac[0] = 1.0;

    // Utilize exp(-zeta*(p+1)**2) = exp(-zeta*p**2)*exp(-2*zeta*p)*exp(-zeta)
% for d in range(1, order + 1):
    ffac[${d}] = ffac[${d - 1}]*v*v*f;
    v *= f;
% endfor

    // Compute filtered solution
    for (int uidx = 0; uidx < ${nupts}; uidx++)
    {
        for (int vidx = 0; vidx < ${nvars}; vidx++)
        {
            fpdtype_t tmp = 0.0;

            // Group terms by basis order
        % for d in range(order + 1):
            tmp += ffac[${d}]*(${' + '.join(f'vdm[uidx][{k}]*umodes[{k}][vidx]'
                                              for k, dd in enumerate(ubdegs) if dd == d)});
        % endfor

            uf[uidx][vidx] = tmp;
        }
    }
</%pyfr:macro>

<%pyfr:macro name='apply_filter_single' params='up, f, rho, rhoY, p, e'>

    fpdtype_t u[${nvars}];

    // Start accumulation
    % for vidx in range(nvars):
        u[${vidx}] = up[0][${vidx}];
    % endfor

    // Apply filter to local value
    fpdtype_t v = 1.0, v2 = 1.0;
    for (int pidx = 1; pidx < ${order+1}; pidx++)
    {
        // Utilize exp(-zeta*(p+1)**2) = exp(-zeta*p**2)*exp(-2*zeta*p)*exp(-zeta)
        v2 *= v*v*f;
        v *= f;

        % for vidx in range(nvars):
        u[${vidx}] += v2*up[pidx][${vidx}];
        % endfor
    }
    rho = u[0];
    ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};
    ${pyfr.expand('get_min_rhoY', 'u', 'q', 'rhoY')};
    p = q[0];
    ${pyfr.expand('compute_entropy', 'u', 'q', 'e')};

</%pyfr:macro>

<%pyfr:kernel name='entropyfilter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='inout fpdtype_t[${str(nfaces)}]'
              vdm='in broadcast fpdtype_t[${str(nefpts)}][${str(nupts)}]'
              invvdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'
              sensor='in fpdtype_t[2]'
              zeta='out fpdtype_t'>

    fpdtype_t rhomin, rhoYmin, pmin, emin;

    // Compute minimum entropy from current and adjacent elements
    fpdtype_t entmin = ${fpdtype_max};
    for (int fidx = 0; fidx < ${nfaces}; fidx++) entmin = fmin(entmin, entmin_int[fidx]);

    // Check if solution is within bounds
    ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'rhoYmin', 'pmin', 'emin')};

    // Filter if out of bounds
    if (rhomin < ${d_min} || rhoYmin < 0.0 || pmin < ${p_min} || emin < entmin - ${e_tol})
    {
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
        fpdtype_t f = 1.0;
        fpdtype_t f_low, f_high, fnew;

        fpdtype_t rho, rhoY, p, e;

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
            ${pyfr.expand('apply_filter_single', 'up', 'f', 'rho', 'rhoY', 'p', 'e')};

            // Update f if constraints aren't satisfied
            if (rho < ${d_min} || rhoY < 0.0 || p < ${p_min} || e < entmin - ${e_tol})
            {
                // Set root-finding interval
                f_high = f;
                f_low = 0.0;

                // Compute brackets
                ${pyfr.expand('apply_filter_single', 'up', 'f_low', 'rho', 'rhoYmin', 'p', 'e')};

                for (int iter = 0; iter < ${niters} && f_high - f_low > ${f_tol}; iter++)
                {
                    // Compute new guess using bisection
                    fnew = 0.5*(f_low + f_high);

                    // Compute filtered state
                    ${pyfr.expand('apply_filter_single', 'up', 'fnew', 'rho', 'rhoY', 'p', 'e')};

                    // Update brackets
                    if (rho < ${d_min} || rhoY < 0.0 || p < ${p_min} || e < entmin - ${e_tol})
                        f_high = fnew;
                    else
                        f_low = fnew;
                }

                // Set current minimum f as the bounds-preserving value
                f = f_low;
            }
        }
        // Filter full solution with bounds-preserving f value
        ${pyfr.expand('apply_filter_full', 'umodes', 'vdm', 'u', 'f')};

        // Calculate minimum entropy from filtered solution
        ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'rhoYmin', 'pmin', 'emin')};
    }

    // Set new minimum entropy within element for next stage
% for fidx in range(nfaces):
    entmin_int[${fidx}] = emin;
% endfor

</%pyfr:kernel>
