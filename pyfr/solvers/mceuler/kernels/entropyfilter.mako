<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.entropy'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.intestar'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\

<%pyfr:macro name='get_minima' params='u, m0, rhomin, tot_rhoYmin, rhoYmin, intemin, emin, Xmin, entmin'>

    rhomin = ${fpdtype_max};
    tot_rhoYmin = ${fpdtype_max};
    % for n in range(ns):
    rhoYmin[${n}] = ${fpdtype_max};
    % endfor
    intemin = ${fpdtype_max};
    emin = ${fpdtype_max};
    Xmin = ${fpdtype_max};

    fpdtype_t ui[${nvars}];
    fpdtype_t qi[${nvars + 2}];
    fpdtype_t qhi[${4 + ns}];

    for (int i = 0; i < ${nupts}; i++)
    {

    % for j in range(nvars):
        ui[${j}] = u[i][${j}];
    % endfor

        // Compute thermodynamic properties
        ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

        fpdtype_t e;
        ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};

        fpdtype_t intestar;
        ${pyfr.expand('compute_intestar', 'ui', 'qi', 'qhi', 'intestar')};

        rhomin = fmin(rhomin, qi[${rhoix}]);
        % for n in range(ns):
          rhoYmin[${n}] = fmin(rhoYmin[${n}], ui[${n}]);
          tot_rhoYmin = fmin(tot_rhoYmin, ui[${n}]);
        % endfor
        intemin = fmin(intemin, intestar);
        emin = fmin(emin, e);
        Xmin = fmin(Xmin, qi[${rhoix}]*(e - entmin));
    }

    % if not fpts_in_upts:
    for (int fidx = 0; fidx < ${nfpts}; fidx++)
    {
        % for vidx in range(nvars):
        ui[${vidx}] = ${pyfr.dot('m0[fidx][{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
        % endfor

        // Compute thermodynamic properties
        ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

        fpdtype_t e;
        ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};

        fpdtype_t intestar;
        ${pyfr.expand('compute_intestar', 'ui', 'qi', 'qhi', 'intestar')};

        rhomin = fmin(rhomin, qi[${rhoix}]);
        % for n in range(ns):
          rhoYmin[${n}] = fmin(rhoYmin[${n}], ui[${n}]);
          tot_rhoYmin = fmin(tot_rhoYmin, ui[${n}]);
        % endfor
        intemin = fmin(intemin, intestar);
        emin = fmin(emin, e);
        Xmin = fmin(Xmin, qi[${rhoix}]*(e - entmin));
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

<%pyfr:macro name='apply_filter_single' params='up, f, rho, rhoY, inte, e, X, entmin'>

    fpdtype_t ui[${nvars}];
    fpdtype_t qi[${nvars + 2}];
    fpdtype_t qhi[${4 + ns}];

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
    ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};
    rho = qi[${rhoix}];
    rhoY = ${fpdtype_max};
    % for n in range(ns):
      rhoY = fmin(rhoY, ui[${n}]);
    % endfor
    X = rho*(e - entmin);

</%pyfr:macro>

<%pyfr:kernel name='entropyfilter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='inout fpdtype_t[${str(nfaces)}]'
              vdm='in broadcast fpdtype_t[${str(nefpts)}][${str(nupts)}]'
              invvdm='in broadcast fpdtype_t[${str(nupts)}][${str(nupts)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'>

    fpdtype_t rhomin, tot_rhoYmin, rhoYmin[${ns}], intemin, emin, Xmin;

    // Compute minimum entropy from current and adjacent elements
    fpdtype_t entmin = ${fpdtype_max};
    for (int fidx = 0; fidx < ${nfaces}; fidx++) entmin = fmin(entmin, entmin_int[fidx]);

    // Check if solution is within bounds
    ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'emin', 'Xmin', 'entmin')};

    // Filter if out of bounds
    if (rhomin < ${d_min} || tot_rhoYmin < 0.0 || intemin < ${inte_min} || Xmin < ${-e_tol})
    {
        % if linearise:

        // Compute mean state
        fpdtype_t uavg[${nvars}], rhoYavg[${ns}], intestaravg, eavg;
        % for vidx in range(nvars):
        uavg[${vidx}] = ${' + '.join(f'{jx}*u[{j}][{vidx}]'
                                     for j, jx in enumerate(meanwts) if jx != 0)};
        % endfor

        fpdtype_t qavg[${nvars + 2}];
        fpdtype_t qhavg[${4 + ns}];
        ${pyfr.expand('stateFrom-cons', 'uavg', 'qavg', 'qhavg')};
        ${pyfr.expand('compute_intestar', 'uavg', 'qavg', 'qhavg', 'eavg')};
        ${pyfr.expand('compute_entropy', 'uavg', 'qavg', 'eavg')};

        fpdtype_t Xavg = qavg[${rhoix}]*(eavg - entmin);


        // Apply density, species, internal energy, and entropy limiting sequentially
        // Density positivity
        if (rhomin < ${d_min}){
            fpdtype_t alpha = (rhomin - ${d_min})/(rhomin - qavg[${rhoix}]);
            alpha = fmin(fmax(alpha, 0.0), 1.0);

            % for uidx, sidx in pyfr.ndrange(nupts, ns):
            u[${uidx}][${sidx}] += alpha*(uavg[${sidx}] - u[${uidx}][${sidx}]);
            % endfor

            // Get new updated values
            ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'emin', 'Xmin', 'entmin')};
        }

        // Species mass >= 0
        % for n in range(ns):
        {
            if (rhoYmin[${n}] < 0.0){
                fpdtype_t rYmin = rhoYmin[${n}];
                fpdtype_t alpha = (rYmin - 0.0)/(rYmin - uavg[${n}]);
                alpha = fmin(fmax(alpha, 0.0), 1.0);

                % for uidx in range(nupts):
                u[${uidx}][${n}] += alpha*(uavg[${n}] - u[${uidx}][${n}]);
                % endfor

                // Get new updated values
                ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'emin', 'Xmin', 'entmin')};
            }
        }
        % endfor

        // Shifted internal energy positivity
        if (intemin < ${inte_min}){
            fpdtype_t alpha = (intemin - ${inte_min})/(intemin - intestaravg);
            alpha = fmin(fmax(alpha, 0.0), 1.0);

            % for uidx, vidx in pyfr.ndrange(nupts, nvars):
            u[${uidx}][${vidx}] += alpha*(uavg[${vidx}] - u[${uidx}][${vidx}]);
            % endfor

            // Get new updated values
            ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'emin', 'Xmin', 'entmin')};
        }

        // Apply minimum entropy principle X = r*(s - s0)
        if (Xmin < ${-e_tol}){
            // Filter to X = 0 (not tolerance)
            fpdtype_t alpha = (Xmin - (0.0))/(Xmin - Xavg);
            alpha = fmin(fmax(alpha, 0.0), 1.0);

            % for uidx, vidx in pyfr.ndrange(nupts, nvars):
            u[${uidx}][${vidx}] += alpha*(uavg[${vidx}] - u[${uidx}][${vidx}]);
            % endfor

            // Get new updated values
            ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'emin', 'Xmin', 'entmin')};
        }

        // Non-linear filtering
        % else:
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

        fpdtype_t rho, rhoY, inte, e, X;

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
            ${pyfr.expand('apply_filter_single', 'up', 'f', 'rho', 'rhoY', 'inte', 'e', 'X', 'entmin')};

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
                    ${pyfr.expand('apply_filter_single', 'up', 'fnew', 'rho', 'rhoY', 'inte', 'e', 'X', 'entmin')};

                    // Update brackets
                    if (rho < ${d_min} || rhoY < 0.0 || inte < ${inte_min} || X < ${-e_tol})
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
        ${pyfr.expand('get_minima', 'u', 'm0', 'rhomin', 'tot_rhoYmin', 'rhoYmin', 'intemin', 'emin', 'Xmin', 'entmin')};
        %endif
    }

    // Set new minimum entropy within element for next stage
% for fidx in range(nfaces):
    entmin_int[${fidx}] = emin;
% endfor
</%pyfr:kernel>
