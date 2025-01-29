<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## <%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

## % if bccfluxstate:
## <%include file='pyfr.solvers.navstokes.kernels.bcs.${bccfluxstate}'/>
## % endif

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ul='in view fpdtype_t[${str(nfpts)}][${str(nvars)}]'>

    printf("ELEMENT\n");
    for(int i=0; i < ${nupts}; i++)
    {
      printf("%03.0f %03.0f %03.0f %03.0f \n", u[i][0], u[i][1], u[i][2], u[i][3]);
    }
    printf("FACE\n");
    % for fpt in range(nfacefpts):
      printf("%03.0f %03.0f %03.0f %03.0f \n", ul[${facefpts[fpt]}][0], ul[${facefpts[fpt]}][1], ul[${facefpts[fpt]}][2], ul[${facefpts[fpt]}][3]);
    % endfor

</%pyfr:kernel>
