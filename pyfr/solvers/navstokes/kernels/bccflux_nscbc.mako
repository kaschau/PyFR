<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>
<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

<%def name="nidx(comp,phys)">
  <% return comp*ndims + phys %>
</%def>\
<% from math import sqrt %>\

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='inout view fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              nl='in fpdtype_t[${str(nfacefpts)}][${str(ndims)}]'
              smats_u='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              jacs='in fpdtype_t[${str(nfacefpts)}]'>

printf("\n*************ELEMENT************\n");

## Step 1: Compute transformed visc flux, pressure, velocity at solution points
fpdtype_t tFi[${nupts}][${ndims}][${nvars}] = {{{0}}};
fpdtype_t dtFidEupts[${nupts}][${ndims}][${nvars}] = {{{0}}};
fpdtype_t pupts[${nupts}];
fpdtype_t vupts[${nupts}][${ndims}];
for (int uidx = 0; uidx < ${nupts}; uidx++)
{
  fpdtype_t ul[${nvars}];
  % for var in range(nvars):
    ul[${var}] = u[uidx][${var}];
  % endfor
  fpdtype_t fi[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  ${pyfr.expand('inviscid_flux', 'ul', 'fi', 'p', 'v')};
  pupts[uidx] = p;
  % for dim in range(ndims):
    vupts[uidx][${dim}] = v[${dim}];
  % endfor

  % for var in range(nvars):
    % for comp in range(ndims):
      % for phys in range(ndims):
        tFi[uidx][${comp}][${var}] += smats_u[uidx][${nidx(comp,phys)}]*fi[${phys}][${var}];
      % endfor
    % endfor
  % endfor
}

% for upt in range(nupts):
  % for upt2 in range(nupts):
    % for var in range(nvars):
      % for comp in range(ndims):
          dtFidEupts[${upt}][${comp}][${var}] += tFi[${upt2}][${comp}][${var}]*${m1[upt,comp,upt2]};
      % endfor
    % endfor
  % endfor
% endfor
## Check
## % for var in range(nvars):
## % for upt in range(nupts):
##   ## printf("upt ${upt} tFi_${var} = %.1f %.1fe \n", tFi[${upt}][0][${var}], tFi[${upt}][1][${var}]);
##   printf("upt ${upt} dtFidEupts_${var} = %e %e \n", dtFidEupts[${upt}][0][${var}], dtFidEupts[${upt}][1][${var}]);
## % endfor
## % endfor

## Iterate over the flux points on our face
% for f, fpt_idx in enumerate(facefpts):
{
  printf("\n Flux point %d\n", ${fpt_idx});

  ## Check
  ## % for var in range(nvars):
  ##   printf("uf ${var} = %f \n", uf[${fpt_idx}][${var}]);
  ## % endfor

  ## Step 1a: Compute physical flux point values
  fpdtype_t f_f[${ndims}][${nvars}];
  fpdtype_t f_f_n[${nvars}] = {0};
  fpdtype_t p_f, v_f[${ndims}];
  fpdtype_t u_f[${nvars}];
  % for var in range(nvars):
    u_f[${var}] = uf[${fpt_idx}][${var}];
  % endfor
  ${pyfr.expand('inviscid_flux', 'u_f', 'f_f', 'p_f', 'v_f')};

  % for upt in range(nupts):
  % for comp in range(ndims):
  % for var in range(nvars):
    f_f_n[${var}] += tFi[${upt}][${comp}][${var}]*${m2[f,comp,upt]};
  % endfor
  % endfor
  % endfor

  ## Check
  ## % for var in range(nvars):
  ##   printf("f_f ${var} = %f %f \n", f_f[0][${var}], f_f[1][${var}]);
  ##   printf("f_f_n ${var} = %f \n", f_f_n[${var}]);
  ## % endfor
  ##   printf("p_f = %f  vf = %f %f \n", p_f, v_f[0], v_f[1]);

  ## Step 2: Compute derivative in transformed space of tflux, pressure, velocity at flux point
  ## Also compute gradient of smats for geometric source term
  fpdtype_t dtFidE[${ndims}][${nvars}] = {{0}};
  fpdtype_t dsmatsdE[${ndims}][${ndims}] = {{0}};
  fpdtype_t dvdE[${ndims}][${ndims}] = {{0}};
  fpdtype_t dpdE[${ndims}] = {0};
  fpdtype_t drhodE[${ndims}] = {0};
  % for upt in range(nupts):
    % for comp in range(ndims):
      % for var in range(nvars):
        ## dtFidE[${comp}][${var}] += tFi[${upt}][${comp}][${var}]*${m12[f, comp, upt]};
        dtFidE[${comp}][${var}] += dtFidEupts[${upt}][${comp}][${var}]*${m0[f, upt]};
      % endfor
      % for phys in range(ndims):
        dvdE[${phys}][${comp}] += vupts[${upt}][${phys}]*${m12[f, comp, upt]};
        dsmatsdE[${comp}][${phys}] += smats_u[${upt}][${nidx(comp,phys)}]*${m12[f, comp, upt]};
      % endfor
      dpdE[${comp}] += pupts[${upt}]*${m12[f, comp, upt]};
      drhodE[${comp}] += u[${upt}][0]*${m12[f, comp, upt]};
    % endfor
  % endfor

  ## Check
  ## % for var in range(nvars):
  ##   printf("dtFidE_${var} = %.1f %.1f \n", dtFidE[0][${var}],dtFidE[1][${var}]);
  ## % endfor
  ## printf("dudE = %f %f\n", dvdE[0][0], dvdE[0][1]);
  ## printf("dvdE = %f %f\n", dvdE[1][0], dvdE[1][1]);
  ## printf("dExdE = %f %f\n", dsmats[0][0], dsmats[0][1]);
  ## printf("dEydE = %f %f\n", dsmats[1][0], dsmats[1][1]);
  ## printf("dpdE = %f %f\n", dpdE[0], dpdE[1]);
  ## printf("drhodE = %f %f\n", drhodE[0], drhodE[1]);

  ## Step 3: Realign derivatives to a face-normal orientation where \Xi is normal to face
  fpdtype_t bnorm[${ndims}] = {${','.join([str(i) for i in bnorm_facefpts[f,:]])}};
  fpdtype_t dtFidE_n[${ndims}][${nvars}];
  fpdtype_t dsmatsdE_n[${ndims}][${nvars}];
  % for var in range(nvars):
  {
    fpdtype_t dF_temp[${ndims}] = {${','.join([f'dtFidE[{i}][{var}]' for i in range(ndims)])}};
    fpdtype_t dF_n_temp[${ndims}];
    ${pyfr.expand('transform_to', 'bnorm', 'dF_temp', 'dF_n_temp', off=0)};
    % for comp in range(ndims):
      dtFidE_n[${comp}][${var}] = dF_n_temp[${comp}];
    % endfor
  }
  % endfor
  fpdtype_t dvdE_n[${ndims}][${ndims}];
  % for phys in range(ndims):
  {
    fpdtype_t dv_temp[${ndims}] = {${','.join([f'dvdE[{phys}][{comp}]' for comp in range(ndims)])}};
    fpdtype_t dv_n_temp[${ndims}];
    ${pyfr.expand('transform_to', 'bnorm', 'dv_temp', 'dv_n_temp', off=0)};
    fpdtype_t dsmats_temp[${ndims}] = {${','.join([f'dsmatsdE[{comp}][{phys}]' for comp in range(ndims)])}};
    fpdtype_t dsmats_n_temp[${ndims}];
    ${pyfr.expand('transform_to', 'bnorm', 'dsmats_temp', 'dsmats_n_temp', off=0)};
    % for comp in range(ndims):
      dvdE_n[${phys}][${comp}] = dv_n_temp[${comp}];
      dsmatsdE_n[${comp}][${phys}] = dsmats_n_temp[${comp}];
    % endfor
  }
  % endfor
  fpdtype_t dpdE_n[${ndims}];
  ${pyfr.expand('transform_to', 'bnorm', 'dpdE', 'dpdE_n', off=0)};
  fpdtype_t drhodE_n[${ndims}];
  ${pyfr.expand('transform_to', 'bnorm', 'drhodE', 'drhodE_n', off=0)};

  ## Check
  ## % for var in range(nvars):
  ##   printf("dtFidE_n_${var} = %f %f \n", dtFidE_n[0][${var}],dtFidE_n[1][${var}]);
  ## % endfor
  ## printf("dudE_n = %f %f\n", dvdE_n[0][0], dvdE_n[0][1]);
  ## printf("dvdE_n = %f %f\n", dvdE_n[1][0], dvdE_n[1][1]);
  ## printf("dpdE_n = %f %f\n", dpdE_n[0], dpdE_n[1]);
  ## printf("drhodE_n = %f %f\n", drhodE_n[0], drhodE_n[1]);

  ## Step 4: Compute the characteristic wave amplitudes, L
  fpdtype_t nl_f[${ndims}] = {${", ".join([f'nl[{f}][{i}]' for i in range(ndims)])}};
  fpdtype_t mag_nl = sqrt(${pyfr.dot('nl_f[{i}]', i=ndims)});
  fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl_f[{i}]', i=ndims)};

  fpdtype_t uhat = (1.0/jacs[${f}])*(${pyfr.dot('nl_f[{i}]','v_f[{i}]', i=ndims)});
  fpdtype_t c = sqrt(${c['gamma']}*p_f/u_f[0]);
  fpdtype_t csq = c*c;
  fpdtype_t chat = c*mag_nl/jacs[${f}];
  fpdtype_t rcpcsq = 1.0/csq;

  fpdtype_t L[${nvars}];

  % if ndims == 2:
    L[0] = uhat*norm_nl[0]*(drhodE_n[0] - dpdE_n[0]*rcpcsq);
    L[1] = uhat*norm_nl[1]*(drhodE_n[0] - dpdE_n[0]*rcpcsq);
    L[2] = ${1.0/sqrt(2)}*(uhat + chat)*(norm_nl[0]*dvdE_n[0][0] + norm_nl[1]*dvdE_n[1][0] + 1.0/(u_f[0]*c)*dpdE_n[0]);
    L[3] = ${1.0/sqrt(2)}*(uhat - chat)*(-norm_nl[0]*dvdE_n[0][0] - norm_nl[1]*dvdE_n[1][0] + 1.0/(u_f[0]*c)*dpdE_n[0]);
  % endif

  ## Check
  % for var in range(nvars):
    printf("L_${var} %f\n", L[${var}]);
  % endfor
  printf("nl_f %f %f \n", nl_f[0], nl_f[1]);
  printf("norm_nl %f %f \n", norm_nl[0], norm_nl[1]);
  printf("v_f %f %f \n", v_f[0], v_f[1]);
  printf("mag_nl %f p_f %f rho_f %f\n", mag_nl, p_f, u_f[0]);
  printf("uhat %f chat %f \n", uhat, chat);

  ## Step 5: Replace incoming wave amplitudes
  ## ${pyfr.expand('compute_L', 'L', 'u_f')};
  fpdtype_t Msq = (${pyfr.dot('v_f[{i}]', i=ndims)})*rcpcsq;
  L[3] = jacs[${f}]*${c['sigma']/sqrt(2)}/u_f[0]*(1.0-Msq)*(p_f - ${c['p']});

  ## Check
  % for i in range(nvars):
    printf("L_${i} = %f \n", L[${i}]);
  % endfor

  ## Step 6: Compute d,dFstardE values normal to face
  fpdtype_t d[${nvars}];
  % if ndims == 2:
    d[0] = norm_nl[0]*L[0] + norm_nl[1]*L[1] + u_f[0]/(${sqrt(2)}*c)*(L[2] + L[3]);
    d[1] = norm_nl[0]*${1.0/sqrt(2)}*(L[2] - L[3]);
    d[2] = norm_nl[1]*${1.0/sqrt(2)}*(L[2] - L[3]);
    d[3] = u_f[0]*${1.0/sqrt(2)}*c*(L[2] + L[3]);

    dtFidE_n[0][0] = jacs[${f}]*(d[0]);
    dtFidE_n[0][1] = jacs[${f}]*(v_f[0] * d[0] + u_f[0]*d[1]);
    dtFidE_n[0][2] = jacs[${f}]*(v_f[1] * d[0] + u_f[0]*d[2]);
    dtFidE_n[0][3] = jacs[${f}]*(0.5*(${pyfr.dot('v_f[{i}]', i=ndims)})*d[0] + u_f[1]*d[1] + u_f[2]*d[2] + ${1.0/(c['gamma']-1.0)}*d[3]);
  % endif

  ## Check
  fpdtype_t dtFidE_star[${ndims}][${nvars}];
  % for var in range(nvars):
  {
    fpdtype_t dF_n_temp[${ndims}] = {${','.join([f'dtFidE_n[{i}][{var}]' for i in range(ndims)])}};
    fpdtype_t dF_temp[${ndims}];
    ${pyfr.expand('transform_from', 'bnorm', 'dF_n_temp', 'dF_temp', off=0)};;
    % for comp in range(ndims):
      dtFidE_star[${comp}][${var}] = dF_temp[${comp}];
    % endfor
  }
  % endfor

  ## Step 7: Solve for normal transformed common flux

  % for var in range(nvars):
  {
    int var = ${var};
    ## printf("\n VAR ${var} \n");
    ## we have dudt (~\del \dot ~f) at our flux point
    uf[${fpt_idx}][${var}] = ${'+'.join([f'dtFidE_star[{dim}][{var}]' for dim in range(ndims)])};
    printf("dtFidE_* %.14e \n", ${'+'.join([f'dtFidE_star[{dim}][{var}]' for dim in range(ndims)])});

    ## Add geomteric source term
    ## uf[${fpt_idx}][${var}] += (${pyfr.dot('f_f[{i}][var]','dsmatsdE_n[0][{i}]', i=ndims)});
    ## printf("geom source %f\n",${pyfr.dot('f_f[{i}][var]','dsmatsdE_n[0][{i}]', i=ndims)});

    ## subtract our flux gradient on the face (from interior values)
    uf[${fpt_idx}][${var}] -= (${'+'.join([f'dtFidE[{dim}][{var}]' for dim in range(ndims)])});
    ## printf("dtFidE %f \n", (${'+'.join([f'dtFidE[{dim}][{var}]' for dim in range(ndims)])}));

    ## divide by ~\del \dot g
    uf[${fpt_idx}][${var}] /= ${m11[f]};

    ## compute and add the transformed, normal flux from interior values
    ## fpdtype_t f_f_n = ${pyfr.dot('nl_f[{i}]', 'f_f[{i}][var]', i=ndims)};
    ## printf("f_f_n  %f\n", f_f_n);
    uf[${fpt_idx}][${var}] += f_f_n[${var}];
    ## Check
    printf("f_n  %f\n", uf[${fpt_idx}][${var}]);
  }
  % endfor

}
% endfor

</%pyfr:kernel>
