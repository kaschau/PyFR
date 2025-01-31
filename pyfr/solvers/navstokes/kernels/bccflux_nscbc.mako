<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

## <%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

## % if bccfluxstate:
## <%include file='pyfr.solvers.navstokes.kernels.bcs.${bccfluxstate}'/>
## % endif

<%def name="nidx(comp,phys)">
  <% return comp*ndims + phys %>
</%def>\

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='inout view fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              nl='in fpdtype_t[${str(nfacefpts)}][${str(ndims)}]'
              smats='in fpdtype_t[${str(nfacefpts)}][${str(ndims*ndims)}]'>

    ## printf("\n*************ELEMENT************\n");
    ## for(int i=0; i < ${nupts}; i++)
    ## {
    ##   printf("%.1f %.1f %.1f %.1f \n", u[i][0], u[i][1], u[i][2], u[i][3]);
    ## }
    ## fpdtype_t dfdE[${str(nvars)}];
    ## fpdtype_t dfdN[${str(nvars)}];
    ## printf("FACE\n");
    ## % for i,fpt in enumerate(facefpts):
    ##   printf("ul = %.1f %.1f %.1f %.1f \n", ul[${fpt}][0], ul[${fpt}][1], ul[${fpt}][2], ul[${fpt}][3]);
    ##   printf("nl = %.1f %.1f\n", nl[${i}][0], nl[${i}][1]);

    ##   % for n in range(nvars):
    ##     dfdE[${n}] = 0.0;
    ##     dfdN[${n}] = 0.0;
    ##   % endfor
    ##   % for upt in range(nupts):
    ##     dfdE[0] += u[${upt}][0]*${m12[fpt,0,upt]};
    ##     dfdE[1] += u[${upt}][1]*${m12[fpt,0,upt]};
    ##     dfdE[2] += u[${upt}][2]*${m12[fpt,0,upt]};
    ##     dfdE[3] += u[${upt}][3]*${m12[fpt,0,upt]};

    ##     dfdN[0] += u[${upt}][0]*${m12[fpt,1,upt]};
    ##     dfdN[1] += u[${upt}][1]*${m12[fpt,1,upt]};
    ##     dfdN[2] += u[${upt}][2]*${m12[fpt,1,upt]};
    ##     dfdN[3] += u[${upt}][3]*${m12[fpt,1,upt]};
    ##   % endfor
    ##   printf("dudE %.1f %.1f %.1f %.1f \n", dfdE[0], dfdE[1], dfdE[2], dfdE[3]);
    ##   printf("dedN %.1f %.1f %.1f %.1f \n", dfdN[0], dfdN[1], dfdN[2], dfdN[3]);

    ## % for comp in range(ndims):
    ## % for phys in range(ndims):
    ##   printf("d${comp}/d${phys} fp = %.1f \n", smats[${i}][${nidx(comp,phys)}]);
    ## % endfor
    ## % endfor

    ## % endfor

    fpdtype_t fl[${ndims}][${nvars}];
    fpdtype_t ul[${nvars}];
    fpdtype_t v_t[${ndims}];
    fpdtype_t p, v[${ndims}];
    fpdtype_t eta_x[${ndims}];
    fpdtype_t c, c_t;

    ## Iterate over the flux points on our face
    % for f,fpt in enumerate(facefpts):

    ## Compute flux based on interior data
    % for i in range(nvars):
      ul[${i}] = uf[${fpt}][${i}];
    % endfor
    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'p', 'v')};

    ## Physical velocity to transformed velocity
    % for i in range(ndims):
      v_t[${i}] = 0;
      % for j in range(ndims):
        <% n = nidx(i,j) %>
        v_t[${i}] += smats[${f}][${n}]*v[${i}];
      % endfor
    % endfor

    ## Compute the transformation vector normal to our transformed face
    % for i in range(ndims):
      eta_x[${i}] = 0;
      % for j in range(ndims):
        <% n = nidx(i,j) %>
        eta_x[${i}] += smats[${f}][${n}]*${bnorm_fpts[i]};
      % endfor
    % endfor

    ## Compute speed of sound normal to the face
    c = sqrt(${c['gamma']}*p/ul[0]);
    c_t = c*sqrt(${pyfr.dot('eta_x[{i}]', i=ndims)});

    % endfor
</%pyfr:kernel>
