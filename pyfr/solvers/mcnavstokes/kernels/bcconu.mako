<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.mixture_state'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.rhoe-from-rhoTY'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${trans}'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bcconu' ndim='1'
              ulin='in view fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'
              nlin='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nlin[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nlin[{i}]', i=ndims)};

    ${pyfr.expand('bc_ldg_state', 'ulin', 'norm_nl', 'ulout')};
</%pyfr:kernel>
