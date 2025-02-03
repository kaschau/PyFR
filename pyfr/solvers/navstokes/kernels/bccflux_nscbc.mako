<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>
<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

<%def name="nidx(comp,phys)">
  <% return comp*ndims + phys %>
</%def>\

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='inout view fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              nl='in fpdtype_t[${str(nfacefpts)}][${str(ndims)}]'
              smats_f='in fpdtype_t[${str(nfacefpts)}][${str(ndims*ndims)}]'
              smats_u='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'>

printf("\n*************ELEMENT************\n");

## Compute transformed visc flux at solution points
fpdtype_t Fi[${nupts}][${ndims}][${nvars}] = {{{0}}};
fpdtype_t pu[${nupts}];
for (int uidx = 0; uidx < ${nupts}; uidx++)
{
  fpdtype_t ul[${nvars}];
  % for vidx in range(nvars):
    ul[${vidx}] = u[uidx][${vidx}];
  % endfor
  fpdtype_t fi_temp[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  pu[uidx] = p;
  ${pyfr.expand('inviscid_flux', 'ul', 'fi_temp', 'p', 'v')};

  % for vidx in range(nvars):
    % for comp in range(ndims):
      % for phys in range(ndims):
        Fi[uidx][${comp}][${vidx}] += smats_u[uidx][${nidx(comp,phys)}]*fi_temp[${phys}][${vidx}];
      % endfor
    % endfor
  % endfor
}

## Iterate over the flux points on our face
% for f, fpt_idx in enumerate(facefpts):
{
  printf("Flux point %d\n", ${fpt_idx});

  ## Step 1: Transform smats to a face-normal orientation where \Xi is normal to face
  fpdtype_t bnorm[${ndims}] = {${','.join([str(i) for i in bnorm_facefpts[f,:]])}};
  fpdtype_t smatsf_t[${ndims}][${ndims}];
  % for phys in range(ndims):
  {
    fpdtype_t smatsf_temp[${ndims}] = {${','.join([f'smats_f[{f}][{nidx(comp,phys)}]' for comp in range(ndims)])}};
    fpdtype_t smatsf_t_temp[${ndims}];
    ${pyfr.expand('transform_to', 'bnorm', 'smatsf_temp', 'smatsf_t_temp', off=0)};
    % for comp in range(ndims):
      smatsf_t[${comp}][${phys}] = smatsf_t_temp[${comp}];
    % endfor
  }
  % endfor

  ## Step 2a: Compute transformed (regular transformed coords) transformed flux derivatives at flux point based
  fpdtype_t dFdE[${ndims}][${nvars}] = {{0}};
  % for dim in range(ndims):
    % for vidx in range(nvars):
      % for upt in range(nupts):
        dFdE[${dim}][${vidx}] += Fi[${upt}][${dim}][${vidx}]*${m12[fpt_idx, dim, upt]};
      % endfor
    % endfor
  % endfor

  ## CHECK
  ## % for i in range(ndims):
  ##   % for j in range(nvars):
  ##     printf("dFdE[${i}][${j}] = %.1f\n", dFdE[${i}][${j}]);
  ##   % endfor
  ## % endfor

  ## Step 2b: Convert these flux derivatives into into the face normal transformed coordinates
  fpdtype_t dFdE_t[${ndims}][${nvars}];
  % for vidx in range(nvars):
  {
    fpdtype_t dF_temp[${ndims}] = {${','.join([f'dFdE[{i}][{vidx}]' for i in range(ndims)])}};
    fpdtype_t dF_t_temp[${ndims}];
    ${pyfr.expand('transform_to', 'bnorm', 'dF_temp', 'dF_t_temp', off=0)};
    % for i in range(ndims):
      dFdE_t[${i}][${vidx}] = dF_t_temp[${i}];
    % endfor
  }
  % endfor

  ## CHECK
  ## % for i in range(ndims):
  ##   % for j in range(nvars):
  ##     printf("dFdE_t[${i}][${j}] = %.1f\n", dFdE_t[${i}][${j}]);
  ##   % endfor
  ## % endfor

  ## Step 3: Compute initial guess of



  ## ## Physical velocity to transformed velocity
  ## fpdtype_t v_t[${ndims}];
  ## % for i in range(ndims):
  ##   v_t[${i}] = 0;
  ##   % for j in range(ndims):
  ##     <% n = nidx(i,j) %>
  ##     v_t[${i}] += smats[${f}][${n}]*v[${i}];
  ##   % endfor
  ## % endfor

  ## ## Compute the transformation vector normal to our transformed face
  ## fpdtype_t eta_x[${ndims}];
  ## % for i in range(ndims):
  ##   eta_x[${i}] = 0;
  ##   % for j in range(ndims):
  ##     <% n = nidx(i,j) %>
  ##     eta_x[${i}] += smats[${f}][${n}]*${bnorm_facefpts[i]};
  ##   % endfor
  ## % endfor

  ## ## Compute speed of sound normal to the face
  ## fpdtype_t c, c_t;
  ## c = sqrt(${c['gamma']}*p/ul[0]);
  ## c_t = c*sqrt(${pyfr.dot('eta_x[{i}]', i=ndims)});

}
% endfor

</%pyfr:kernel>
