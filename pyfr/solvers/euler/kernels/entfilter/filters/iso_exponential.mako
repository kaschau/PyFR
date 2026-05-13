<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entropy'/>

// Isotropic exponential modal filter:  sigma_k = f^(d_k^2)
//
// f = 1 is the identity; f -> 0 collapses everything onto mode 0 (the cell
// mean).  The rolling-power trick exploits exp(-zeta*(p+1)^2) =
// exp(-zeta*p^2)*exp(-2*zeta*p)*exp(-zeta) to avoid per-mode pow() calls.

<%pyfr:macro name='iso_exp_apply_full' params='umodes, vdm, uf, f'>
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

<%pyfr:macro name='iso_exp_apply_single' params='up, f, d, p, e'>
    // Start accumulation
    fpdtype_t ui[${nvars}];
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

    ${pyfr.expand('compute_entropy', 'ui', 'd', 'p', 'e' )};
</%pyfr:macro>
