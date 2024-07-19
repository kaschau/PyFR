<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-prims'/>
<%include file='pyfr.solvers.mceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mceuler.kernels.bcs.${bctype}'/>

##

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};
<% ns = c['ns'] %>

    // Compute left thermodynamic quantities
    fpdtype_t ql[${nvars + 1}];
    fpdtype_t qhl[${3 + ns}];
    ${pyfr.expand('stateFrom-cons', 'ul', 'ql', 'qhl')};

    // Compute the right BC state
    fpdtype_t ur[${nvars}];
    fpdtype_t qr[${nvars + 1}];
    fpdtype_t qhr[${3 + ns}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'ql', 'qhl', 'norm_nl', 'ur', 'qr', 'qhr')};

    // Perform the Riemann solve
    fpdtype_t fn[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'ql', 'qr', 'qhl', 'qhr', 'norm_nl', 'fn')};

    // Scale and write out the common normal fluxes
% for i in range(nvars):
    ul[${i}] = mag_nl*fn[${i}];
% endfor
</%pyfr:kernel>
