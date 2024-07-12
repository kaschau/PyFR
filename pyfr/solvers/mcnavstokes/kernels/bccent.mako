<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.${bctype}'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.entropy'/>

<%pyfr:kernel name='bccent' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              entmin_lhs='out view reduce(min) fpdtype_t'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    // Compute left thermodynamic quantities
    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('stateFrom-cons', 'ul', 'ql', 'qhl')};

    // Compute the right BC state
    fpdtype_t ur[${nvars}];
    fpdtype_t qr[${nvars+1}];
    fpdtype_t qhr[${3+ns}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'ql', 'qhl', 'norm_nl', 'ur', 'qr', 'qhr')};

    // Compute entropy for boundary state
    fpdtype_t p, d, entmin_rhs;
    ${pyfr.expand('compute_entropy', 'ur', 'qr', 'entmin_rhs')};

    // Compute face minima (reduce with atomics)
    entmin_lhs = entmin_rhs;
</%pyfr:kernel>
