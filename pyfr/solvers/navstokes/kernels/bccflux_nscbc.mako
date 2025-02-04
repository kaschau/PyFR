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
              smats_f='in fpdtype_t[${str(nfacefpts)}][${str(ndims*ndims)}]'
              smats_u='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              jacs='in fpdtype_t[${str(nfacefpts)}]'>

printf("\n*************ELEMENT************\n");

## Step 1: Compute transformed visc flux, pressure, velocity at solution points
fpdtype_t tFi[${nupts}][${ndims}][${nvars}] = {{{0}}};
fpdtype_t pupts[${nupts}];
fpdtype_t vupts[${nupts}][${ndims}];
for (int uidx = 0; uidx < ${nupts}; uidx++)
{
  fpdtype_t ul[${nvars}];
  % for vidx in range(nvars):
    ul[${vidx}] = u[uidx][${vidx}];
  % endfor
  fpdtype_t fi[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  ${pyfr.expand('inviscid_flux', 'ul', 'fi', 'p', 'v')};
  pupts[uidx] = p;
  % for dim in range(ndims):
    vupts[uidx][${dim}] = v[${dim}];
  % endfor

  % for vidx in range(nvars):
    % for comp in range(ndims):
      % for phys in range(ndims):
        tFi[uidx][${comp}][${vidx}] += smats_u[uidx][${nidx(comp,phys)}]*fi[${phys}][${vidx}];
      % endfor
    % endfor
  % endfor
}
## Check
## % for vidx in range(nvars):
## % for upt in range(nupts):
##   printf("upt ${upt} tFi_${vidx} = %.1f %.1fe \n", tFi[${upt}][0][${vidx}], tFi[${upt}][1][${vidx}]);
## % endfor
## % endfor

## Iterate over the flux points on our face
% for f, fpt_idx in enumerate(facefpts):
{
  printf("Flux point %d\n", ${fpt_idx});

  ## Step 1a: Compute physical flux point values
  fpdtype_t f_f[${ndims}][${nvars}];
  fpdtype_t p_f, v_f[${ndims}];
  fpdtype_t u_f[${nvars}];
  % for vidx in range(nvars):
    u_f[${vidx}] = uf[${fpt_idx}][${vidx}];
  % endfor
  ${pyfr.expand('inviscid_flux', 'u_f', 'f_f', 'p_f', 'v_f')};

  ## Step 2: Compute derivatives in transformed space of tflux, pressure, velocity at flux point
  fpdtype_t dtFidE[${ndims}][${nvars}] = {{0}};
  fpdtype_t dvdE[${ndims}][${ndims}] = {{0}};
  fpdtype_t dpdE[${ndims}] = {0};
  fpdtype_t drhodE[${ndims}] = {0};
  % for comp in range(ndims):
    % for upt in range(nupts):
      % for vidx in range(nvars):
        dtFidE[${comp}][${vidx}] += tFi[${upt}][${comp}][${vidx}]*${m12[f, comp, upt]};
      % endfor
      % for phys in range(ndims):
        dvdE[${phys}][${comp}] += vupts[${upt}][${phys}]*${m12[f, comp, upt]};
      % endfor
      dpdE[${comp}] += pupts[${upt}]*${m12[f, comp, upt]};
      drhodE[${comp}] += u[${upt}][0]*${m12[f, comp, upt]};
    % endfor
  % endfor

  ## Check
  ## % for vidx in range(nvars):
  ##   printf("dtFidE_${vidx} = %.1f %.1f \n", dtFidE[0][${vidx}],dtFidE[1][${vidx}]);
  ## % endfor
  ## printf("dudE = %.1f %.1f\n", dvdE[0][0], dvdE[0][1]);
  ## printf("dvdE = %.1f %.1f\n", dvdE[1][0], dvdE[1][1]);
  ## printf("dpdE = %.1f %.1f\n", dpdE[0], dpdE[1]);
  ## printf("drhodE = %.1f %.1f\n", drhodE[0], drhodE[1]);


  ## Step 3: Realign derivatives to a face-normal orientation where \Xi is normal to face
  fpdtype_t bnorm[${ndims}] = {${','.join([str(i) for i in bnorm_facefpts[f,:]])}};
  fpdtype_t dtFidE_n[${ndims}][${nvars}];
  % for vidx in range(nvars):
  {
    fpdtype_t dF_temp[${ndims}] = {${','.join([f'dtFidE[{i}][{vidx}]' for i in range(ndims)])}};
    fpdtype_t dF_n_temp[${ndims}];
    ${pyfr.expand('transform_to', 'bnorm', 'dF_temp', 'dF_n_temp', off=0)};
    % for comp in range(ndims):
      dtFidE_n[${comp}][${vidx}] = dF_n_temp[${comp}];
    % endfor
  }
  % endfor
  fpdtype_t dvdE_n[${ndims}][${ndims}];
  % for phys in range(ndims):
  {
    fpdtype_t dv_temp[${ndims}] = {${','.join([f'dvdE[{phys}][{comp}]' for comp in range(ndims)])}};
    fpdtype_t dv_n_temp[${ndims}];
    ${pyfr.expand('transform_to', 'bnorm', 'dv_temp', 'dv_n_temp', off=0)};
    % for comp in range(ndims):
      dvdE_n[${phys}][${comp}] = dv_n_temp[${comp}];
    % endfor
  }
  % endfor
  fpdtype_t dpdE_n[${ndims}];
  ${pyfr.expand('transform_to', 'bnorm', 'dpdE', 'dpdE_n', off=0)};
  fpdtype_t drhodE_n[${ndims}];
  ${pyfr.expand('transform_to', 'bnorm', 'drhodE', 'drhodE_n', off=0)};

  ## Check
  ## % for vidx in range(nvars):
  ##   printf("dtFidE_n_${vidx} = %.1f %.1f \n", dtFidE_n[0][${vidx}],dtFidE_n[1][${vidx}]);
  ## % endfor
  ## printf("dudE_n = %.1f %.1f\n", dvdE_n[0][0], dvdE_n[0][1]);
  ## printf("dvdE_n = %.1f %.1f\n", dvdE_n[1][0], dvdE_n[1][1]);
  ## printf("dpdE_n = %.1f %.1f\n", dpdE_n[0], dpdE_n[1]);
  ## printf("drhodE_n = %.1f %.1f\n", drhodE_n[0], drhodE_n[1]);

  ## Step 4: Compute the characteristic wave amplitudes, L
  fpdtype_t nl_f[${ndims}] = {${", ".join([f'nl[{f}][{i}]' for i in range(ndims)])}};
  fpdtype_t mag_nl = sqrt(${pyfr.dot('nl_f[{i}]', i=ndims)});
  fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl_f[{i}]', i=ndims)};

  fpdtype_t uhat = ${pyfr.dot('nl_f[{i}]','v_f[{i}]', i=ndims)};
  fpdtype_t c = sqrt(${c['gamma']}*p_f/u_f[0]);
  fpdtype_t csq = c*c;
  fpdtype_t chat = c*mag_nl;
  fpdtype_t rcpcsq = 1.0/csq;

  fpdtype_t L[${nvars}];

  % if ndims == 2:
    L[0] = uhat*norm_nl[0]*(drhodE_n[0]- dpdE[0]*rcpcsq);
    L[1] = uhat*norm_nl[1]*(drhodE_n[0] - dpdE[0]*rcpcsq);
    L[2] = ${1.0/sqrt(2)}*(uhat + chat)*(norm_nl[0]*dvdE_n[0][0] + norm_nl[1]*dvdE_n[1][0] + 1.0/(u_f[0]*c)*dpdE[0]);
    L[3] = ${1.0/sqrt(2)}*(uhat - chat)*(-norm_nl[0]*dvdE_n[0][0] - norm_nl[1]*dvdE_n[1][0] + 1.0/(u_f[0]*c)*dpdE[0]);
  % endif

  ## Check
  ## printf("norm_nl %f %f \n", norm_nl[0], norm_nl[1]);
  ## printf("v_f %f %f \n", v_f[0], v_f[1]);
  ## printf("mag_nl %f p_f %f rho_f %f\n", mag_nl, p_f, u_f[0]);
  ## printf("uhat %f chat %f \n", uhat, chat);

  ## Step 5: Replace incoming wave amplitudes
  ## ${pyfr.expand('compute_L', 'L', 'u_f')};
  fpdtype_t Msq = ${pyfr.dot('v_f[{i}]', i=ndims)}*rcpcsq;
  L[3] = jacs[${f}]*${c['sigma']/sqrt(2)/c['Lx']}/u_f[0]*(1.0-Msq)*(p_f - ${c['p_inf']});

  ## Check
  ## % for i in range(nvars):
  ##   printf("L_${i} = %f \n", L[${i}]);
  ## % endfor

  ## Step 6: Compute dFstardE values normal to face
  % if ndims == 2:
    dtFidE_n[0][0] = norm_nl[0]*L[0] + norm_nl[1]*L[1] + u_f[0]/(${sqrt(2)}*c)*(L[2] + L[3]);
    dtFidE_n[0][1] = norm_nl[0]*${1.0/sqrt(2)}*(L[2] - L[3]);
    dtFidE_n[0][2] = norm_nl[1]*${1.0/sqrt(2)}*(L[2] - L[3]);
    dtFidE_n[0][3] = u_f[0]*${1.0/sqrt(2)}*(L[2] + L[3]);
  % endif

  ## Step 7: Solve for normal transformed common flux

  % for vidx in range(nvars):
  {
    ## we have dudt (~\del \dot ~f) at our flux point
    uf[${fpt_idx}][${vidx}] = ${'+'.join([f'dtFidE_n[{dim}][{vidx}]' for dim in range(ndims)])};

    ## subtract our flux gradient on the face (from interior values)
    uf[${fpt_idx}][${vidx}] -= (${'+'.join([f'dtFidE[{dim}][{vidx}]' for dim in range(ndims)])});

    ## divide by ~\del \dot g
    uf[${fpt_idx}][${vidx}] /= ${m11[f]};

    ## compute and add the transformed, normal flux from interior values
    int idx = ${vidx};
    fpdtype_t f_f_n = ${pyfr.dot('nl_f[{i}]', 'f_f[{i}][idx]', i=ndims)};
    uf[${fpt_idx}][${vidx}] += f_f_n;

    ## Check
    printf("f_n %f\n", uf[${fpt_idx}][${vidx}]);
  }
  % endfor

  ## ## Transform smats to a face-normal orientation where \Xi is normal to face
  ## fpdtype_t bnorm[${ndims}] = {${','.join([str(i) for i in bnorm_facefpts[f,:]])}};
  ## fpdtype_t smatsf_t[${ndims}][${ndims}];
  ## % for phys in range(ndims):
  ## {
  ##   fpdtype_t smatsf_temp[${ndims}] = {${','.join([f'smats_f[{f}][{nidx(comp,phys)}]' for comp in range(ndims)])}};
  ##   fpdtype_t smatsf_t_temp[${ndims}];
  ##   ${pyfr.expand('transform_to', 'bnorm', 'smatsf_temp', 'smatsf_t_temp', off=0)};
  ##   % for comp in range(ndims):
  ##     smatsf_t[${comp}][${phys}] = smatsf_t_temp[${comp}];
  ##   % endfor
  ## }
  ## % endfor

  ## ## Step 2a: Compute transformed (regular transformed coords) transformed flux derivatives at flux point based
  ## fpdtype_t dFdE[${ndims}][${nvars}] = {{0}};
  ## % for dim in range(ndims):
  ##   % for vidx in range(nvars):
  ##     % for upt in range(nupts):
  ##       dFdE[${dim}][${vidx}] += Fi[${upt}][${dim}][${vidx}]*${m12[f, dim, upt]};
  ##     % endfor
  ##   % endfor
  ## % endfor

  ## ## CHECK
  ## ## % for i in range(ndims):
  ## ##   % for j in range(nvars):
  ## ##     printf("dFdE[${i}][${j}] = %.1f\n", dFdE[${i}][${j}]);
  ## ##   % endfor
  ## ## % endfor

  ## ## Step 2b: Convert these flux derivatives into into the face normal transformed coordinates
  ## fpdtype_t dFdE_t[${ndims}][${nvars}];
  ## % for vidx in range(nvars):
  ## {
  ##   fpdtype_t dF_temp[${ndims}] = {${','.join([f'dFdE[{i}][{vidx}]' for i in range(ndims)])}};
  ##   fpdtype_t dF_t_temp[${ndims}];
  ##   ${pyfr.expand('transform_to', 'bnorm', 'dF_temp', 'dF_t_temp', off=0)};
  ##   % for i in range(ndims):
  ##     dFdE_t[${i}][${vidx}] = dF_t_temp[${i}];
  ##   % endfor
  ## }
  ## % endfor

  ## ## CHECK
  ## ## % for i in range(ndims):
  ## ##   % for j in range(nvars):
  ## ##     printf("dFdE_t[${i}][${j}] = %.1f\n", dFdE_t[${i}][${j}]);
  ## ##   % endfor
  ## ## % endfor

  ## ## Step 3: Compute initial guess of



  ## ## ## Physical velocity to transformed velocity
  ## ## fpdtype_t v_t[${ndims}];
  ## ## % for i in range(ndims):
  ## ##   v_t[${i}] = 0;
  ## ##   % for j in range(ndims):
  ## ##     <% n = nidx(i,j) %>
  ## ##     v_t[${i}] += smats[${f}][${n}]*v[${i}];
  ## ##   % endfor
  ## ## % endfor

  ## ## ## Compute the transformation vector normal to our transformed face
  ## ## fpdtype_t eta_x[${ndims}];
  ## ## % for i in range(ndims):
  ## ##   eta_x[${i}] = 0;
  ## ##   % for j in range(ndims):
  ## ##     <% n = nidx(i,j) %>
  ## ##     eta_x[${i}] += smats[${f}][${n}]*${bnorm_facefpts[i]};
  ## ##   % endfor
  ## ## % endfor

  ## ## ## Compute speed of sound normal to the face
  ## ## fpdtype_t c, c_t;
  ## ## c = sqrt(${c['gamma']}*p/ul[0]);
  ## ## c_t = c*sqrt(${pyfr.dot('eta_x[{i}]', i=ndims)});

}
% endfor

</%pyfr:kernel>
