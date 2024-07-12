<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.entropy'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='get_minmax_Y' params='q, Ymin, Ymax'>

    Ymin = ${fpdtype_max};
    Ymax = -${fpdtype_max};

    % for n in range(ns):
    Ymin = fmin(Ymin, q[${Yix + n}]);
    Ymax = fmax(Ymax, q[${Yix + n}]);
    % endfor
</%pyfr:macro>

<%pyfr:macro name='get_minima' params='u, Ymin, Ymax, pmin, emin'>
    fpdtype_t ui[${nvars}];

    Ymin = ${fpdtype_max};
    Ymax = -${fpdtype_max};
    pmin = ${fpdtype_max};
    emin = ${fpdtype_max};

    for (int i = 0; i < ${nupts}; i++)
    {
    % for j in range(nvars):
        ui[${j}] = u[i][${j}];
    % endfor

        // Compute thermodynamic properties
        fpdtype_t qi[${nvars+1}];
        fpdtype_t qhi[${3+ns}];
        ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

        fpdtype_t e;
        ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};
        % for n in range(ns):
        Ymin = fmin(Ymin, qi[${Yix + n}]);
        Ymax = fmax(Ymax, qi[${Yix + n}]);
        % endfor
        pmin = fmin(pmin, qi[0]);
        emin = fmin(emin, e);
    }
</%pyfr:macro>

<%pyfr:macro name='apply_filter_full' params='umodes, vdm, uf, f'>
    // Precompute filter factors per basis degree
    fpdtype_t ffac[${order + 1}];
    fpdtype_t v = ffac[0] = 1.0;

    // Utilize exp(-zeta*p**2) = pow(f, p**2)
% for d in range(1, order + 1):
    v *= f;
    ffac[${d}] = v*v;
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

<%pyfr:macro name='apply_filter_single' params='u, q, qh, f'>
    // Apply filter to local value
    fpdtype_t v = 1.0;
    for (int pidx = 1; pidx < ${order+1}; pidx++)
    {
        // Utilize exp(-zeta*p**2) = pow(f, p**2)
        v *= f;

        % for vidx in range(nvars):
        u[${vidx}] += v*v*up[pidx][${vidx}];
        % endfor
    }
    ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};
</%pyfr:macro>

<%pyfr:kernel name='entropyfilter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='inout fpdtype_t[${str(nfaces)}]'
              vdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              invvdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'>

    fpdtype_t Ymin, Ymax, pmin, emin;

    // Compute minimum entropy from current and adjacent elements
    fpdtype_t entmin = ${fpdtype_max};
    for (int fidx = 0; fidx < ${nfaces}; fidx++) entmin = fmin(entmin, entmin_int[fidx]);

    // Check if solution is within bounds
    ${pyfr.expand('get_minima', 'u', 'Ymin', 'Ymax', 'pmin', 'emin')};

    // Filter if out of bounds
    if (Ymin < ${Y_min} || Ymax > ${Y_max} || pmin < ${p_min} || emin < entmin - ${e_tol})
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

        fpdtype_t Y_min, Y_max, p, e;
        fpdtype_t Y_min_low, Y_max_low, p_low, e_low;
        fpdtype_t Y_min_high, Y_max_high, p_high, e_high;

        // Compute f on a rolling basis per solution point
        fpdtype_t up[${order+1}][${nvars}];
        for (int uidx = 0; uidx < ${nupts}; uidx++)
        {
            // Group nodal contributions by common filter factor
        % for pidx, vidx in pyfr.ndrange(order+1, nvars):
            up[${pidx}][${vidx}] = (${' + '.join(f'vdm[uidx][{k}]*umodes[{k}][{vidx}]'
                                                   for k, dd in enumerate(ubdegs) if dd == pidx)});
        % endfor

            // Compute constraints with current minimum f value
            fpdtype_t ui[${nvars}];
            // Compute thermodynamic properties
            fpdtype_t qi[${nvars+1}];
            fpdtype_t qhi[${3+ns}];
            ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

            // Start accumulation
            % for vidx in range(nvars):
                ui[${vidx}] = up[0][${vidx}];
            % endfor

            ${pyfr.expand('apply_filter_single', 'ui', 'qi', 'qhi', 'f')};
            ${pyfr.expand('get_minmax_Y', 'qi', 'Ymin', 'Ymax')};
            p = qi[0];
            ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};


            // Update f if constraints aren't satisfied
            if (Ymin < ${Y_min} || Ymax > ${Y_max} || p < ${p_min} || e < entmin - ${e_tol})
            {
                // Set root-finding interval
                f_high = f;
                f_low = 0.0;

                // Compute brackets
                Y_min_high = Ymin; Y_max_high = Ymax;
                p_high = p; e_high = e;
                ${pyfr.expand('apply_filter_single', 'ui', 'qi', 'qhi', 'f_low')};

                ${pyfr.expand('get_minmax_Y', 'ui', 'Y_min_low', 'Y_max_low')};
                p_low = qi[0];
                ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e_low')};

                // Regularize constraints to be around zero
                Y_min_low -= ${Y_min}; Y_min_high -= ${Y_min};
                Y_max_low -= ${Y_max}; Y_max_high -= ${Y_min};
                p_low -= ${p_min}; p_high -= ${p_min};
                e_low -= entmin - ${e_tol}; e_high -= entmin - ${e_tol};

                // Iterate filter strength with bisection algorithm
                for (int iter = 0; iter < ${niters} && f_high - f_low > ${f_tol}; iter++)
                {
                    // Compute new guess using bisection
                    fnew = 0.5*(f_low + f_high);

                    // Compute filtered state
                    ${pyfr.expand('apply_filter_single', 'ui', 'qi', 'qhi', 'fnew')};

                    ${pyfr.expand('get_minmax_Y', 'ui', 'Y_min', 'Y_max')};
                    p = qi[0];
                    ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};

                    // Update brackets
                    if (Ymin < ${Y_min} || Ymax > ${Y_max} || p < ${p_min} || e < entmin - ${e_tol})
                    {
                        f_high = fnew;
                        Y_min_high = Ymin - ${Y_min};
                        Y_max_high = Ymax - ${Y_max};
                        p_high = p - ${p_min};
                        e_high = e - (entmin - ${e_tol});
                    }
                    else
                    {
                        f_low = fnew;
                        Y_min_low = Ymin - ${Y_min};
                        Y_max_low = Ymax - ${Y_max};
                        p_low = p - ${p_min};
                        e_low = e - (entmin - ${e_tol});
                    }
                }

                // Set current minimum f as the bounds-preserving value
                f = f_low;
            }
        }

        // Filter full solution with bounds-preserving f value
        ${pyfr.expand('apply_filter_full', 'umodes', 'vdm', 'u', 'f')};

        // Calculate minimum entropy from filtered solution
        ${pyfr.expand('get_minima', 'u', 'Ymin', 'Ymax', 'pmin', 'emin')};
    }

    // Set new minimum entropy within element for next stage
% for fidx in range(nfaces):
    entmin_int[${fidx}] = emin;
% endfor
</%pyfr:kernel>
