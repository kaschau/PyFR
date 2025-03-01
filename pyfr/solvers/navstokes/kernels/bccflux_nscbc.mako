<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>
<%include file='pyfr.solvers.navstokes.kernels.flux'/>
<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

<%def name="nidx(comp,phys)">
  <% return comp*ndims + phys %>
</%def>\
<%def name="nfidx(dim, fpt)">
  <% return dim*nfpts + fpt %>
</%def>\
<%def name="nuidx(dim, upt)">
  <% return dim*nupts + upt %>
</%def>\
<%from math import sqrt %>\

<%pyfr:macro name='PU.dE' params='dE, N'>
fpdtype_t nx = norm_nl[0];
fpdtype_t ny = norm_nl[1];

fpdtype_t u_ = v[0];
fpdtype_t v_ = v[1];
fpdtype_t gmo = ${c['gamma'] - 1.0};

fpdtype_t k = 0.5*(${pyfr.dot('v[{i}]', i=ndims)});

fpdtype_t dE0 = dE[0];
fpdtype_t dE1 = dE[1];
fpdtype_t dE2 = dE[2];
fpdtype_t dE4 = dE[3];

N[0] = dE0 + gmo*(-dE4 - dE0*k + dE1*u_ + dE2*v_)*invcsq;
N[1] = (dE2*nx-dE1*ny+dE0*ny*u_-dE0*nx*v_)*invrho;
N[2] = ${1.0/sqrt(2)}*(dE4*gmo + c*dE1*nx + c*dE2*ny - dE1*gmo*u_ - dE2*gmo*v_ + dE0*(gmo*k - c*nx*u_ - c*ny*v_))*invc*invrho;
N[3] = ${1.0/sqrt(2)}*(dE4*gmo - c*dE1*nx - c*dE2*ny - dE1*gmo*u_ - dE2*gmo*v_ + dE0*(gmo*k + c*nx*u_ + c*ny*v_))*invc*invrho;
</%pyfr:macro>


<%pyfr:macro name='PUinv.N' params='N, dE'>
fpdtype_t nx = norm_nl[0];
fpdtype_t ny = norm_nl[1];

fpdtype_t u_ = v[0];
fpdtype_t v_ = v[1];
fpdtype_t gmo = ${c['gamma'] - 1.0};

fpdtype_t k = 0.5*(${pyfr.dot('v[{i}]', i=ndims)});

fpdtype_t N0 = N[0];
fpdtype_t N1 = N[1];
fpdtype_t Np = N[2];
fpdtype_t Nm = N[3];

dE[0] = N0 + ${1.0/sqrt(2)}*(Nm+Np)*rho*invc;
dE[1] = -N1*ny*rho + N0*u_ + ${1.0/sqrt(2)}*rho*invc*(Nm*(-c*nx+u_) + Np*(c*nx+u_));
dE[2] = N1*nx*rho + ${1.0/sqrt(2)}*rho*(-Nm*ny + Np*ny) + N0*v_ + ${1.0/sqrt(2)}*rho*invc*v_*(Nm+Np);
dE[3] = k*N0 + ${1.0/sqrt(2)}*(c*(Nm+Np)*rho/gmo + k*(Nm+Np)*rho*invc)-0.5*rho*(${sqrt(2)}*(Nm-Np)*(nx*u_+ny*v_) + 2.0*N1*(ny*u_-nx*v_));
</%pyfr:macro>

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u_upts='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              u_fpts='inout view fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              gradu_upts='in view fpdtype_t[${str(ndims*nupts)}][${str(nvars)}]'
              gradu_fpts='in view fpdtype_t[${str(ndims*nfpts)}][${str(nvars)}]'
              normnl_ffpt='in fpdtype_t[${str(nfacefpts)}][${str(ndims)}]'
              smats_upts='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              jacs_ffpt='in fpdtype_t[${str(nfacefpts)}]'>

## printf("\n*************ELEMENT************\n");

## Step 1: Compute transformed and physical flux at solution points
fpdtype_t f_upts[${nupts}][${ndims}][${nvars}] = {{{0}}};
fpdtype_t tF_upts[${nupts}][${ndims}][${nvars}] = {{{0}}};
% for upt in range(nupts):
{
  fpdtype_t u[${nvars}];
  fpdtype_t gradu[${ndims}][${nvars}];
  % for var in range(nvars):
    u[${var}] = u_upts[${upt}][${var}];
    % for dim in range(ndims):
      gradu[${dim}][${var}] = gradu_upts[${nuidx(dim,upt)}][${var}];
    % endfor
  % endfor
  fpdtype_t f[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  ${pyfr.expand('inviscid_flux', 'u', 'f', 'p', 'v')};

  ## TODO: cache blocking can eliminate the gradu_upts
  ## so need to address to incorporate viscous fluxes
  ## ${pyfr.expand('viscous_flux_add', 'u', 'gradu', f')};


  % for var in range(nvars):
    % for phys in range(ndims):
      % for comp in range(ndims):
        tF_upts[${upt}][${comp}][${var}] += smats_upts[${upt}][${nidx(comp,phys)}]*f[${phys}][${var}];
      % endfor
      f_upts[${upt}][${phys}][${var}] = f[${phys}][${var}];
    % endfor
  % endfor
}
% endfor

## Check
## % for upt in range(nupts):
##   % for var in range(nvars):
##     printf("upt ${upt} tFi_${var} = %.14e %.14e \n", tFi[${upt}][0][${var}], tFi[${upt}][1][${var}]);
##   % endfor
## % endfor

## % for upt in range(nupts):
## % for var in range(nvars):
##   printf("upt ${upt} ${var} = %.2f \n", u_upts[${upt}][${var}]);
## % endfor
## % endfor

## % for upt in range(nupts):
## % for var in range(nvars):
##   printf("grad_upt ${upt} ${var} = %.2f %.2f \n", gradu_upts[${nuidx(0,upt)}][${var}], gradu_upts[${nuidx(1,upt)}][${var}]);
## % endfor
## % endfor

## Iterate over the flux points on our face
% for f, fpt_idx in enumerate(facefpts):
{
  ## printf("\n Flux point %d\n", ${fpt_idx});

  ## Get face normals at flux point
  fpdtype_t bnorm[${ndims}] = {${','.join([str(i) for i in bnorm_facefpts[f,:]])}};
  fpdtype_t norm_nl[${ndims}] = {${", ".join([f'normnl_ffpt[{f}][{i}]' for i in range(ndims)])}};
  ${pyfr.expand('set_normal', 'bnorm', 'norm_nl')};
  fpdtype_t jac = jacs_ffpt[${f}];

  ## Step 1: Compute transformed flux and metrics relative to face normal
  ## transformed orientation
  fpdtype_t smats_upts_T[${nupts}][${ndims*ndims}]= {{0}};
  fpdtype_t tF_T[${nupts}][${ndims}][${nvars}] = {{{0}}};
  % for upt in range(nupts):
  {
    % for phys in range(ndims):
    {
      fpdtype_t smats_T_temp[${ndims}];
      fpdtype_t smats_temp[${ndims}] = {${','.join([f'smats_upts[{upt}][{nidx(comp,phys)}]' for comp in range(ndims)])}};
      ${pyfr.expand('transform_to', 'bnorm', 'smats_temp', 'smats_T_temp', off=0)};
      % for comp in range(ndims):
        smats_upts_T[${upt}][${nidx(comp,phys)}] = smats_T_temp[${comp}];
      % endfor
    }
    % endfor

    % for var in range(nvars):
      % for comp in range(ndims):
        % for phys in range(ndims):
          tF_T[${upt}][${comp}][${var}] += smats_upts_T[${upt}][${nidx(comp,phys)}]*f_upts[${upt}][${phys}][${var}];
        % endfor
      % endfor
    % endfor
  }
  % endfor

  ## Step 1a: Compute physical flux at our flux point
  fpdtype_t ul[${nvars}];
  fpdtype_t gradul[${ndims}][${nvars}];
  % for var in range(nvars):
    ul[${var}] = u_fpts[${fpt_idx}][${var}];
    % for dim in range(ndims):
      gradul[${dim}][${var}] = gradu_fpts[${nfidx(dim,nfpts)}][${var}];
    % endfor
  % endfor
  fpdtype_t fl[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'p', 'v')};
  ${pyfr.expand('viscous_flux_add', 'ul', 'gradul', 'fl')};

  ## Check
  ## % for var in range(nvars):
  ##   printf("grad_fpt ${fpt_idx} ${var} = %.2f %.2f \n",gradu_fpts[${nfidx(0,fpt_idx)}][${var}], gradu_fpts[${nfidx(1,fpt_idx)}][${var}]);
  ## % endfor

  ## Check
  ## % for var in range(nvars):
  ##   printf("ul ${var} = %f\n", ul[${var}]);
  ##   printf("fl ${var} = %f %f \n", fl[0][${var}], fl[1][${var}]);
  ## % endfor
  ##   printf("p = %f  v = %f %f \n", p, v[0], v[1]);

  ## Compute normal transformed flux
  fpdtype_t fl_n[${nvars}] = {0};
  % for upt in range(nupts):
  % for comp in range(ndims):
  % for var in range(nvars):
    fl_n[${var}] += tF_upts[${upt}][${comp}][${var}]*${m2[f,comp,upt]};
  % endfor
  % endfor
  % endfor

  ## Check
  ## % for var in range(nvars):
  ##   printf("fl ${var} = %f %f \n", fl[0][${var}], fl[1][${var}]);
  ##   printf("fl_n ${var} = %f \n", fl_n[${var}]);
  ## % endfor

  ## Step 3: Compute derivative in normal transformed space of tflux, at flux point
  ## Also compute gradient of smats for geometric source term
  fpdtype_t dtFdE_T[${ndims}][${nvars}] = {{0}};
  fpdtype_t dsmatsdE[${ndims}][${ndims}] = {{0}};
  % for upt in range(nupts):
  {
    fpdtype_t m12_temp[${ndims}] = {${','.join([str(m12[f,i,upt]) for i in range(ndims)])}};
    fpdtype_t m12_T[${ndims}];
    ${pyfr.expand('transform_to', 'bnorm', 'm12_temp', 'm12_T', off=0)};
    % for comp in range(ndims):
      % for var in range(nvars):
        dtFdE_T[${comp}][${var}] += tF_T[${upt}][${comp}][${var}]*m12_T[${comp}];
      % endfor
      % for phys in range(ndims):
        dsmatsdE[${comp}][${phys}] += smats_upts_T[${upt}][${nidx(comp,phys)}]*m12_T[${comp}];
      % endfor
    % endfor
  }
  % endfor

  ## Check
  ## % for var in range(nvars):
  ##   printf("dtFidE_${var} = %f %f \n", dtFidE[0][${var}], dtFidE[1][${var}]);
  ## % endfor

  ## Step 4: Compute the characteristic wave strengths, N
  fpdtype_t rho = ul[0];
  fpdtype_t invrho = 1.0/rho;
  fpdtype_t c = sqrt(${c['gamma']}*p*invrho);
  fpdtype_t invc = 1.0/c;
  fpdtype_t invcsq = invc*invc;
  fpdtype_t N[${nvars}];
  fpdtype_t S[${nvars}]={0};
  fpdtype_t source[${nvars}];
  {
    fpdtype_t dEdE[${nvars}] = {${','.join([f'dtFdE_T[0][{i}]' for i in range(nvars)])}};
    fpdtype_t dGdN[${nvars}] = {${','.join([f'dtFdE_T[1][{i}]' for i in range(nvars)])}};

    % for var in range(nvars):
    {
      source[${var}] = ${'+'.join([f'fl[{dim}][{var}]*dsmatsdE[0][{dim}]' for dim in range(ndims)])};
      dEdE[${var}] -= source[${var}];
      dGdN[${var}] += source[${var}];
    }
    % endfor

    ${pyfr.expand('PU.dE','dEdE','N')}
    ## HACK
    ## ${pyfr.expand('PU.dE','dGdN','S')}
  }

  ## Step 5: Compute wave amplitudes for unknown waves
  ${pyfr.expand('compute_wave_amp', 'ul', 'p', 'v', 'jac', 'N', 'S', 'norm_nl')};

  ## ## Step 6: Compute dtFdE_T* values normal to face
  fpdtype_t dtFdE_Ts[${nvars}];
  {
    ${pyfr.expand('PUinv.N','N','dtFdE_Ts')};
    % for var in range(nvars):
      dtFdE_Ts[${var}] += source[${var}];
    % endfor
  }

  ## Step 7: Solve for normal transformed common flux
  % for var in range(nvars):
  {
    ## printf("\n VAR ${var} \n");
    ## we have dudt (~\del \dot ~f) at our flux point
    u_fpts[${fpt_idx}][${var}] = dtFdE_Ts[${var}];
    ## printf("dtFidE_* %.14e \n", ${'+'.join([f'dtFidE_star[{dim}][{var}]' for dim in range(ndims)])});

    ## subtract our flux gradient on the face (from interior values)
    u_fpts[${fpt_idx}][${var}] -= dtFdE_T[0][${var}];
    ## printf("dtFidE %f \n", (${'+'.join([f'dtFidE[{dim}][{var}]' for dim in range(ndims)])}));

    ## divide by ~\del \dot g
    u_fpts[${fpt_idx}][${var}] /= ${m11[f]};

    ## add the transformed, normal flux from interior values
    u_fpts[${fpt_idx}][${var}] += fl_n[${var}];

    ## Check
    ## printf("f_n =  %.14e \n", u_fpt[${fpt_idx}][${var}]);
  }
  % endfor
}
% endfor

</%pyfr:kernel>
