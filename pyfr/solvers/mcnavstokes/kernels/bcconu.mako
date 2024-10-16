<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-prims'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${eos}.e_Y_Y_x'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${trans}'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.${bctype}'/>

<% ns = c['ns'] %>

<%pyfr:kernel name='bcconu' ndim='1'
              ulin='in view fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'
              nlin='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nlin[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nlin[{i}]', i=ndims)};

    fpdtype_t qlin[${nvars + 1}];
    fpdtype_t qhlin[${4 + ns}];
    ${pyfr.expand('stateFrom-cons', 'ulin', 'qlin', 'qhlin')};
    fpdtype_t qlout[${nvars + 1}];
    fpdtype_t qhlout[${4 + ns}];
    ${pyfr.expand('bc_ldg_state', 'ulin', 'qlin', 'qhlin', 'norm_nl', 'ulout', 'qlout', 'qhlout')};
</%pyfr:kernel>
