<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.entropy'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.intestar'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-cons'/>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>\

<%pyfr:macro name='get_minima' params='u, m0, rhomin, tot_rhoYmin, rhoYmin, intemin, smin, Xmin, s0'>

    rhomin = ${fpdtype_max};
    tot_rhoYmin = ${fpdtype_max};
    % for n in range(mcf.ns):
    rhoYmin[${n}] = ${fpdtype_max};
    % endfor
    intemin = ${fpdtype_max};
    smin = ${fpdtype_max};
    Xmin = ${fpdtype_max};

    fpdtype_t ui[${nvars}];
    fpdtype_t qi[${nvars + 2}];
    fpdtype_t qhi[${4 + mcf.ns}];

    for (int i = 0; i < ${nupts}; i++)
    {

    % for j in range(nvars):
        ui[${j}] = u[i][${j}];
    % endfor

        // Compute thermodynamic properties
        ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

        fpdtype_t s;
        ${pyfr.expand('compute_entropy', 'ui', 'qi', 's')};

        fpdtype_t intestar;
        ${pyfr.expand('compute_intestar', 'ui', 'qi', 'qhi', 'intestar')};

        rhomin = fmin(rhomin, qi[${rhoix}]);
        % for n in range(mcf.ns):
          rhoYmin[${n}] = fmin(rhoYmin[${n}], ui[${n}]);
          tot_rhoYmin = fmin(tot_rhoYmin, ui[${n}]);
        % endfor
        intemin = fmin(intemin, intestar);
        smin = fmin(smin, s);
        Xmin = fmin(Xmin, qi[${rhoix}]*(s - s0));
    }

    % if not fpts_in_upts:
    for (int fidx = 0; fidx < ${nfpts}; fidx++)
    {
        % for vidx in range(nvars):
        ui[${vidx}] = ${pyfr.dot('m0[fidx][{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
        % endfor

        // Compute thermodynamic properties
        ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

        fpdtype_t s;
        ${pyfr.expand('compute_entropy', 'ui', 'qi', 's')};

        fpdtype_t intestar;
        ${pyfr.expand('compute_intestar', 'ui', 'qi', 'qhi', 'intestar')};

        rhomin = fmin(rhomin, qi[${rhoix}]);
        % for n in range(mcf.ns):
          rhoYmin[${n}] = fmin(rhoYmin[${n}], ui[${n}]);
          tot_rhoYmin = fmin(tot_rhoYmin, ui[${n}]);
        % endfor
        intemin = fmin(intemin, intestar);
        smin = fmin(smin, s);
        Xmin = fmin(Xmin, qi[${rhoix}]*(s - s0));
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

<%pyfr:macro name='apply_filter_single' params='up, f, rho, rhoY, inte, s, X, s0'>

    fpdtype_t ui[${nvars}];
    fpdtype_t qi[${nvars + 2}];
    fpdtype_t qhi[${4 + mcf.ns}];

    // Start accumulation
    % for vidx in range(nvars):
        ui[${vidx}] = up[0][${vidx}];
    % endfor

    // Apply filter to local value
    fpdtype_t v = 1.0, v2 = 1.0;
    for (int pidx = 1; pidx < ${order+1}; pidx++)
    {
        // Utilize exp(-zeta*(p+1)**2) = exp(-zeta*p**2)*exp(-2*zeta*p)*exp(-zeta)
        v2 *= v*v*f;
        v *= f;

        % for vidx in range(nvars):
        ui[${vidx}] += v2*up[pidx][${vidx}];
        % endfor
    }

    ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};
    ${pyfr.expand('compute_intestar', 'ui', 'qi', 'qhi', 'inte')};
    ${pyfr.expand('compute_entropy', 'ui', 'qi', 's')};
    rho = qi[${rhoix}];
    rhoY = ${fpdtype_max};
    % for n in range(mcf.ns):
      rhoY = fmin(rhoY, ui[${n}]);
    % endfor
    X = rho*(s - s0);

</%pyfr:macro>

<%pyfr:kernel name='entropyfilter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='inout fpdtype_t[${str(nfaces)}]'
              ef_filter='out fpdtype_t[1]'
              vdm='in broadcast fpdtype_t[${str(nefpts)}][${str(nupts)}]'
              invvdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'
              mean_wts='in fpdtype_t[${str(nupts)}]'>

    fpdtype_t rhomin, tot_rhoYmin, rhoYmin[${mcf.ns}], intemin, smin, Xmin;
    fpdtype_t f = 1.0;

    // Compute minimum entropy from current and adjacent elements
    fpdtype_t s0 = ${fpdtype_max};
    for (int fidx = 0; fidx < ${nfaces}; fidx++) s0 = fmin(s0, entmin_int[fidx]);

    // Check if solution is within bounds
    ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'smin', 'Xmin', 's0')};

    // Compute mean state using per-element weights
    fpdtype_t uavg[${nvars}];
    % for vidx in range(nvars):
    uavg[${vidx}] = ${pyfr.dot('mean_wts[{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
    % endfor

    // Filter if out of bounds
    % if not linearise:
    if (rhomin < ${d_min} || tot_rhoYmin < 0.0 || intemin < ${inte_min} || Xmin < ${-e_tol})
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
        fpdtype_t f_low, f_high, fnew;

        fpdtype_t rho, rhoY, inte, s, X;

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
            ${pyfr.expand('apply_filter_single', 'up', 'f', 'rho', 'rhoY', 'inte', 's', 'X', 's0')};

            // Update f if constraints aren't satisfied
            if (rho < ${d_min} || rhoY < 0.0 || inte < ${inte_min} || X < ${-e_tol})
            {
                // Set root-finding interval
                f_high = f;
                f_low = 0.0;

                for (int iter = 0; iter < ${niters} && f_high - f_low > ${f_tol}; iter++)
                {
                    // Compute new guess using bisection
                    fnew = 0.5*(f_low + f_high);

                    // Compute filtered state
                    ${pyfr.expand('apply_filter_single', 'up', 'fnew', 'rho', 'rhoY', 'inte', 's', 'X', 's0')};

                    // Update brackets
                    if (rho < ${d_min} || rhoY < 0.0 || inte < ${inte_min} || X < ${-e_tol}){
                        f_high = fnew;
                    }else{
                        f_low = fnew;
                    }
                }

                // Set current minimum f as the bounds-preserving value
                f = f_low;
            }
        }

        // Filter full solution with bounds-preserving f value
        ${pyfr.expand('apply_filter_full', 'umodes', 'vdm', 'u', 'f')};

        // Account for mass defect in reference-space filtering
        fpdtype_t duavg[${nvars}] = {0};
        % for vidx in range(nvars):
        duavg[${vidx}] = uavg[${vidx}] - ${pyfr.dot('mean_wts[{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
        % endfor

        % for uidx, vidx in pyfr.ndrange(nupts, nvars):
        u[${uidx}][${vidx}] += duavg[${vidx}];
        % endfor

        // Calculate minimum entropy from filtered solution
        ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'smin', 'Xmin', 's0')};
    }
    % endif

    // Apply linearised limiting
    if (rhomin < ${d_min} || tot_rhoYmin < 0.0 || intemin < ${inte_min} || Xmin < ${-e_tol})
    {
        fpdtype_t qavg[${nvars + 2}];
        fpdtype_t qhavg[${4 + mcf.ns}];
        fpdtype_t intestaravg, savg;
        ${pyfr.expand('stateFrom-cons', 'uavg', 'qavg', 'qhavg')};
        ${pyfr.expand('compute_intestar', 'uavg', 'qavg', 'qhavg', 'intestaravg')};
        ${pyfr.expand('compute_entropy', 'uavg', 'qavg', 'savg')};

        fpdtype_t Xavg = qavg[${rhoix}]*(savg - s0);
        fpdtype_t alpha;

        // Density positivity
        if (rhomin < ${d_min})
        {
            alpha = (rhomin - ${d_min})/(rhomin - qavg[${rhoix}]);
            alpha = fmin(fmax(alpha, 0.0), 1.0);
            f = fmin(f, 1.0 - alpha);

            % for uidx, sidx in pyfr.ndrange(nupts, mcf.ns):
            u[${uidx}][${sidx}] += alpha*(uavg[${sidx}] - u[${uidx}][${sidx}]);
            % endfor

            ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'smin', 'Xmin', 's0')};
        }

        // Species mass >= 0
        fpdtype_t mod = 0.0;
        % for n in range(mcf.ns):
        {
            if (rhoYmin[${n}] < 0.0)
            {
                mod = 1.0;
                fpdtype_t rYmin = rhoYmin[${n}];
                alpha = (rYmin - 0.0)/(rYmin - uavg[${n}]);
                alpha = fmin(fmax(alpha, 0.0), 1.0);
                f = fmin(f, 1.0 - alpha);

                % for uidx in range(nupts):
                u[${uidx}][${n}] += alpha*(uavg[${n}] - u[${uidx}][${n}]);
                % endfor
            }
        }
        % endfor
        if (mod > 0.0)
        {
            ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'smin', 'Xmin', 's0')};
        }

        // Shifted internal energy positivity
        if (intemin < ${inte_min})
        {
            alpha = (intemin - ${inte_min})/(intemin - intestaravg);
            alpha = fmin(fmax(alpha, 0.0), 1.0);
            f = fmin(f, 1.0 - alpha);

            % for uidx, vidx in pyfr.ndrange(nupts, nvars):
            u[${uidx}][${vidx}] += alpha*(uavg[${vidx}] - u[${uidx}][${vidx}]);
            % endfor

            ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'smin', 'Xmin', 's0')};
        }

        // Apply minimum entropy principle X = r*(s - s0)
        if (Xmin < ${-e_tol})
        {
            alpha = (Xmin - (0.0))/(Xmin - Xavg);
            alpha = fmin(fmax(alpha, 0.0), 1.0);
            f = fmin(f, 1.0 - alpha);

            % for uidx, vidx in pyfr.ndrange(nupts, nvars):
            u[${uidx}][${vidx}] += alpha*(uavg[${vidx}] - u[${uidx}][${vidx}]);
            % endfor

            ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'smin', 'Xmin', 's0')};
        }
    }

    // Export filter strength
    ef_filter[0] = f;

    // Set new minimum entropy within element for next stage
% for fidx in range(nfaces):
    entmin_int[${fidx}] = smin;
% endfor
</%pyfr:kernel>
