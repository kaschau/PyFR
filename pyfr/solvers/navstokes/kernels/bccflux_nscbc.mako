<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## <%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

## % if bccfluxstate:
## <%include file='pyfr.solvers.navstokes.kernels.bcs.${bccfluxstate}'/>
## % endif

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'>

    printf("ELEMENT\n");
    for(int i=0; i < ${nupts}; i++)
    {
      printf("%03.0f %03.0f %03.0f %03.0f \n", u[i][0], u[i][1], u[i][2], u[i][3]);
    }

</%pyfr:kernel>
