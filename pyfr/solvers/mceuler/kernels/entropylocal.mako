<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.entropy'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:kernel name='entropylocal' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              entmin_int='out fpdtype_t[${str(nfaces)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'>
    // Compute minimum entropy across element
    fpdtype_t ui[${nvars}], e;

    fpdtype_t entmin = ${fpdtype_max};
    for (int i = 0; i < ${nupts}; i++)
    {
    % for j in range(nvars):
        ui[${j}] = u[i][${j}];
    % endfor

        // Compute thermodynamic properties
        fpdtype_t qi[${nvars + 2}];
        fpdtype_t qhi[${4 + ns}];
        ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

        ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};

        entmin = fmin(entmin, e);
    }
    % if not fpts_in_upts:
    for (int fidx = 0; fidx < ${nfpts}; fidx++)
    {
        % for vidx in range(nvars):
        ui[${vidx}] = ${pyfr.dot('m0[fidx][{k}]', f'u[{{k}}][{vidx}]', k=nupts)};
        % endfor
        // Compute thermodynamic properties
        ${pyfr.expand('stateFrom-cons', 'ui', 'qi', 'qhi')};

        // Compute specific entropy
        ${pyfr.expand('compute_entropy', 'ui', 'qi', 'e')};
        entmin = fmin(entmin, e);
    }
    % endif

    // Set interface entropy values to minimum
% for i in range(nfaces):
    entmin_int[${i}] = entmin;
% endfor
</%pyfr:kernel>
