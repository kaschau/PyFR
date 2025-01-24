<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

% if bccfluxstate:
<%include file='pyfr.solvers.navstokes.kernels.bcs.${bccfluxstate}'/>
% endif

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              ndivg='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};


    printf("HERE!!!! \n");
    printf("%f %f \n", ndivg[0], ndivg[1]);
    printf("%f %f \n", nl[0], nl[1]);
    ${pyfr.expand('bc_common_flux_state', 'ul', 'gradul', 'gradg', 'norm_nl', 'mag_nl')};
</%pyfr:kernel>
