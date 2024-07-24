<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.entropy'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='get_min_Y' params='q, Ymin'>

    Ymin = ${fpdtype_max};
    % for n in range(ns):
    Ymin = fmin(Ymin, q[${Yix + n}]);
    % endfor

</%pyfr:macro>

<%pyfr:macro name='get_minima' params='u, Ymin, rhomin, pmin, emin'>

    Ymin = ${fpdtype_max};
    rhomin = ${fpdtype_max};
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

        Ymin = fmin(Ymin, Ymintemp);
        rhomin = fmin(rhomin, ui[0]);
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

<%pyfr:macro name='apply_filter_single' params='up, f, Ymin, rho, p, e'>

    fpdtype_t u[${nvars}];

    // Start accumulation
    % for vidx in range(nvars):
        u[${vidx}] = up[0][${vidx}];
    % endfor

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

    fpdtype_t q[${nvars+1}];
    fpdtype_t qh[${3+ns}];
    ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};
    ${pyfr.expand('get_min_Y', 'q', 'Ymin')};
    rho = u[0];
    p = q[0];
    ${pyfr.expand('compute_entropy', 'u', 'q', 'e')};

</%pyfr:macro>

<%pyfr:kernel name='entropyfilter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='inout fpdtype_t[${str(nfaces)}]'
              vdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              invvdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              sensor='in fpdtype_t[1][3]'>

    fpdtype_t Ymin, rhomin, pmin, emin;
    fpdtype_t kxrcf = sensor[0][0];

    // Compute minimum entropy from current and adjacent elements
    fpdtype_t entmin = ${fpdtype_max};
    for (int fidx = 0; fidx < ${nfaces}; fidx++) entmin = fmin(entmin, entmin_int[fidx]);

    // Check if solution is within bounds
    ${pyfr.expand('get_minima', 'u', 'Ymin', 'rhomin', 'pmin', 'emin')};

    // Filter if out of bounds
    if (Ymin < -${Y_tol} || rhomin < ${d_min} || pmin < ${p_min} || (emin < entmin - ${e_tol} && kxrcf >= ${s_switch}))
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
        fpdtype_t f1, f2, f3, f4;

        fpdtype_t Ymin, rho, p, e;
        fpdtype_t Ymin_low, rho_low, p_low, e_low;
        fpdtype_t Ymin_high, rho_high, p_high, e_high;

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
            ${pyfr.expand('apply_filter_single', 'up', 'f', 'Ymin', 'rho', 'p', 'e')};

            // Update f if constraints aren't satisfied
            if (Ymin < -${Y_tol} || rho < ${d_min} || p < ${p_min} || (e < entmin - ${e_tol} && kxrcf >= ${s_switch}))
            {
                // Set root-finding interval
                f_high = f;
                f_low = 0.0;

                // Compute brackets
                Ymin_high = Ymin;
                rho_high = rho;
                p_high = p; e_high = e;

                ${pyfr.expand('apply_filter_single', 'up', 'f_low', 'Ymin_low', 'rho_low', 'p_low', 'e_low')};

                // Regularize constraints to be around zero
                Ymin_low += ${Y_tol}; Ymin_high += ${Y_tol};
                rho_low -= ${d_min}; rho_high -= ${d_min};
                p_low -= ${p_min}; p_high -= ${p_min};
                e_low -= entmin - ${e_tol}; e_high -= entmin - ${e_tol};

                // Iterate filter strength with Illinois algorithm
                for (int iter = 0; iter < ${niters} && f_high - f_low > ${f_tol}; iter++)
                {
                    // Compute new guess for each constraint (catch if root is not bracketed)
                    f1 = (Ymin_high > 0.0) ? f_high : (0.5*f_low*Ymin_high - f_high*Ymin_low)/(0.5*Ymin_high - Ymin_low + ${ill_tol});
                    f2 = (rho_high > 0.0) ? f_high : (0.5*f_low*rho_high - f_high*rho_low)/(0.5*rho_high - rho_low + ${ill_tol});
                    f3 = (p_high > 0.0) ? f_high : (0.5*f_low*p_high - f_high*p_low)/(0.5*p_high - p_low + ${ill_tol});
                    f4 = (e_high > 0.0) ? f_high : (0.5*f_low*e_high - f_high*e_low)/(0.5*e_high - e_low + ${ill_tol});

                    // Compute guess as minima of individual constraints
                    fnew = fmin(f1, fmin(f2, fmin(f3, f4)));

                    // In case of bracketing failure (due to roundoff errors), revert to bisection
                    fnew = ((fnew > f_high) || (fnew < f_low)) ? 0.5*(f_low + f_high) : fnew;

                    // Compute filtered state
                    ${pyfr.expand('apply_filter_single', 'up', 'fnew', 'Ymin', 'rho', 'p', 'e')};

                    // Update brackets
                    if (Ymin < -${Y_tol} || rho < ${d_min} || p < ${p_min} || (e < entmin - ${e_tol} && kxrcf >= ${s_switch}))
                    {
                        f_high = fnew;
                        Ymin_high = Ymin + ${Y_tol};
                        rho_high = rho - ${d_min};
                        p_high = p - ${p_min};
                        e_high = e - (entmin - ${e_tol});
                    }
                    else
                    {
                        f_low = fnew;
                        Ymin_low = Ymin + ${Y_tol};
                        rho_low = rho - ${d_min};
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
        ${pyfr.expand('get_minima', 'u', 'Ymin', 'rhomin', 'pmin', 'emin')};
    }

    // Set new minimum entropy within element for next stage
% for fidx in range(nfaces):
    entmin_int[${fidx}] = emin;
% endfor
</%pyfr:kernel>
