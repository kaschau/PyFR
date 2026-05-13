<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.entropy'/>

<%pyfr:macro name='get_minima' params='u, m0, dmin, pmin, emin'>
    fpdtype_t d, p, e;
    fpdtype_t ui[${nvars}];

    dmin = ${fpdtype_max}; pmin = ${fpdtype_max}; emin = ${fpdtype_max};

    for (int i = 0; i < ${nupts}; i++)
    {
    % for j in range(nvars):
        ui[${j}] = u[i][${j}];
    % endfor

        ${pyfr.expand('compute_entropy', 'ui', 'd', 'p', 'e')};
        dmin = fmin(dmin, d); pmin = fmin(pmin, p); emin = fmin(emin, e);
    }

    % if not fpts_in_upts:
    fpdtype_t uf[${nvars}];
    for (int fidx = 0; fidx < ${nfpts}; fidx++)
    {
        % for vidx in range(nvars):
        uf[${vidx}] = ${pyfr.dot('m0[fidx][{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
        % endfor

        ${pyfr.expand('compute_entropy', 'uf', 'd', 'p', 'e')};
        dmin = fmin(dmin, d); pmin = fmin(pmin, p); emin = fmin(emin, e);
    }
    % endif
</%pyfr:macro>
