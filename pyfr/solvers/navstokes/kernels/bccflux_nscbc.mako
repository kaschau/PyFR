<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## <%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

## % if bccfluxstate:
## <%include file='pyfr.solvers.navstokes.kernels.bcs.${bccfluxstate}'/>
## % endif

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ul='in view fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              nl='in fpdtype_t[${str(nfpts)}][${str(ndims)}]'>

    printf("\n*************ELEMENT************\n");
    for(int i=0; i < ${nupts}; i++)
    {
      printf("%.1f %.1f %.1f %.1f \n", u[i][0], u[i][1], u[i][2], u[i][3]);
    }
    fpdtype_t dfdE[${str(nvars)}];
    fpdtype_t dfdN[${str(nvars)}];
    printf("FACE\n");
    % for fpt in facefpts:
      printf("ul = %.1f %.1f %.1f %.1f \n", ul[${fpt}][0], ul[${fpt}][1], ul[${fpt}][2], ul[${fpt}][3]);
      printf("nl = %.1f %.1f\n", nl[${fpt}][0], nl[${fpt}][1]);

      % for n in range(nvars):
        dfdE[${n}] = 0.0;
        dfdN[${n}] = 0.0;
      % endfor
      % for upt in range(nupts):
        dfdE[0] += u[${upt}][0]*${m12[fpt,0,upt]};
        dfdE[1] += u[${upt}][1]*${m12[fpt,0,upt]};
        dfdE[2] += u[${upt}][2]*${m12[fpt,0,upt]};
        dfdE[3] += u[${upt}][3]*${m12[fpt,0,upt]};

        dfdN[0] += u[${upt}][0]*${m12[fpt,1,upt]};
        dfdN[1] += u[${upt}][1]*${m12[fpt,1,upt]};
        dfdN[2] += u[${upt}][2]*${m12[fpt,1,upt]};
        dfdN[3] += u[${upt}][3]*${m12[fpt,1,upt]};
      % endfor
      printf("dudE %.1f %.1f %.1f %.1f \n", dfdE[0], dfdE[1], dfdE[2], dfdE[3]);
      printf("dedN %.1f %.1f %.1f %.1f \n", dfdN[0], dfdN[1], dfdN[2], dfdN[3]);
    % endfor

</%pyfr:kernel>
