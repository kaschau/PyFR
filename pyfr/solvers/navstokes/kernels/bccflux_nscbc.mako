<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>
<%include file='pyfr.solvers.navstokes.kernels.flux'/>

<%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

## <% check = True %>
<% check = False %>

<% sq2 = 2**0.5 %>
<% invsq2 = 2**-0.5 %>
<% pnd = r"%.14e "*ndims %>\

<%def name="nidx(comp,phys)">
  <% return comp*ndims + phys %>
</%def>\
<%def name="nfidx(dim, fpt)">
  <% return dim*nfpts + fpt %>
</%def>\
<%def name="nuidx(dim, upt)">
  <% return dim*nupts + upt %>
</%def>\

<%pyfr:macro name='WU_dot_dE' params='dE, N, u, p, v'>
  fpdtype_t rho = u[0];
  fpdtype_t invrho = 1.0/rho;
  fpdtype_t c = sqrt(${c['gamma']}*p*invrho);
  fpdtype_t invc = 1.0/c;
  fpdtype_t invcsq = invc*invc;

  fpdtype_t nx = norm_nl[0];
  fpdtype_t ny = norm_nl[1];

  fpdtype_t gmo = ${c['gamma'] - 1.0};

  fpdtype_t k = 0.5*(${pyfr.dot('v[{i}]', i=ndims)});

% if ndims == 2:

  N[0] = dE[0] - gmo*invcsq*(dE[0]*k - dE[1]*v[0] - dE[2]*v[1] + dE[3]);
  N[1] = invrho*(dE[0]*(ny*v[0] - nx*v[1]) - dE[1]*ny + dE[2]*nx);
  N[2] = ${invsq2}*invc*invrho*(gmo*(dE[3] + dE[0]*k - (dE[1]*v[0] + dE[2]*v[1])) + c*(dE[1]*nx + dE[2]*ny - dE[0]*(nx*v[0] + ny*v[1])));
  N[3] = ${invsq2}*invc*invrho*(gmo*(dE[3] + dE[0]*k - (dE[1]*v[0] + dE[2]*v[1])) - c*(dE[1]*nx + dE[2]*ny - dE[0]*(nx*v[0] + ny*v[1])));

% elif ndims == 3:

  fpdtype_t nz = norm_nl[2];

  N[0] = gmo*nx*invcsq*(-dE[4] - dE[0]*k + dE[1]*v[0] + dE[2]*v[1] + dE[3]*v[2]) + invrho*(-dE[3]*ny + dE[2]*nz + dE[0]*(nx*rho - nz*v[1] + ny*v[2]));
  N[1] = gmo*ny*invcsq*(-dE[4] - dE[0]*k + dE[1]*v[0] + dE[2]*v[1] + dE[3]*v[2]) + invrho*( dE[3]*nx - dE[1]*nz + dE[0]*(ny*rho + nz*v[0] - nx*v[2]));
  N[2] = gmo*nz*invcsq*(-dE[4] - dE[0]*k + dE[1]*v[0] + dE[2]*v[1] + dE[3]*v[2]) + invrho*(-dE[2]*nx + dE[1]*ny + dE[0]*(nz*rho - ny*v[0] + nx*v[1]));
  N[3] = invc*invrho*${invsq2}*(gmo*(dE[4] + dE[0]*k - (dE[1]*v[0] + dE[2]*v[1] + dE[3]*v[2])) + c*(dE[1]*nx + dE[2]*ny + dE[3]*nz - dE[0]*(nx*v[0] + ny*v[1] + nz*v[2])));
  N[4] = invc*invrho*${invsq2}*(gmo*(dE[4] + dE[0]*k - (dE[1]*v[0] + dE[2]*v[1] + dE[3]*v[2])) - c*(dE[1]*nx + dE[2]*ny + dE[3]*nz - dE[0]*(nx*v[0] + ny*v[1] + nz*v[2])));

% endif
</%pyfr:macro>


<%pyfr:macro name='WUinv_dot_N' params='N, dE, u, p, v'>
  fpdtype_t rho = u[0];
  fpdtype_t invrho = 1.0/rho;
  fpdtype_t c = sqrt(${c['gamma']}*p*invrho);
  fpdtype_t invc = 1.0/c;

  fpdtype_t nx = norm_nl[0];
  fpdtype_t ny = norm_nl[1];

  fpdtype_t gmo = ${c['gamma'] - 1.0};

  fpdtype_t k = 0.5*(${pyfr.dot('v[{i}]', i=ndims)});

  % if ndims == 2:

  dE[0] = N[0] + ${invsq2}*rho*invc*(N[2] + N[3]);
  dE[1] = -N[1]*ny*rho + N[0]*v[0] + ${invsq2}*rho*invc*(N[2]*(c*nx + v[0]) - N[3]*(c*nx - v[0]));
  dE[2] =  N[1]*nx*rho + N[0]*v[1] + ${invsq2}*rho*invc*(N[2]*(c*ny + v[1]) - N[3]*(c*ny - v[1]));
  dE[3] = k*N[0] - N[1]*rho*(ny*v[0] - nx*v[1]) + ${invsq2}*rho*((N[3] + N[2])*(c/gmo + k*invc) + (N[2] - N[3])*(nx*v[0] + ny*v[1]));

  % elif ndims == 3:
  fpdtype_t nz = norm_nl[2];

  dE[0] = N[0]*nx + N[1]*ny + N[2]*nz + ${invsq2}*invc*rho*(N[3] + N[4]);
  dE[1] = rho*(N[2]*ny - N[1]*nz) + N[0]*nx*v[0] + N[1]*ny*v[0] + N[2]*nz*v[0] + ${invsq2}*invc*rho*(N[3]*(c*nx+v[0]) + N[4]*(-c*nx + v[0]));
  dE[2] = rho*(N[0]*nz - N[2]*nx) + N[0]*nx*v[1] + N[1]*ny*v[1] + N[2]*nz*v[1] + ${invsq2}*invc*rho*(N[3]*(c*ny+v[1]) + N[4]*(-c*ny + v[1]));
  dE[3] = rho*(N[1]*nx - N[0]*ny) + N[0]*nx*v[2] + N[1]*ny*v[2] + N[2]*nz*v[2] + ${invsq2}*invc*rho*(N[3]*(c*nz+v[2]) + N[4]*(-c*nz + v[2]));
  dE[4] = N[0]*(k*nx + nz*rho*v[1] - ny*rho*v[2]) +
          N[1]*(k*ny - nz*rho*v[0] + nx*rho*v[2]) +
          N[2]*(k*nz + ny*rho*v[0] - nx*rho*v[1]) +
        rho*invc*${invsq2}/gmo*(N[3]*(c*c + gmo*k + c*gmo*(nx*v[0] + ny*v[1] + nz*v[2]))
                                  + N[4]*(c*c + gmo*k - c*gmo*(nx*v[0] + ny*v[1] + nz*v[2])));

  % endif
</%pyfr:macro>

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u_upts='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              u_fpts='inout view fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              gradu_upts='in view fpdtype_t[${str(ndims*nupts)}][${str(nvars)}]'
              gradu_fpts='in view fpdtype_t[${str(ndims*nfpts)}][${str(nvars)}]'
              normnl_ffpt='in fpdtype_t[${str(nfacefpts)}][${str(ndims)}]'
              smats_upts='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              jacs_ffpt='in fpdtype_t[${str(nfacefpts)}]'>

% if check:
printf("\n*************ELEMENT************\n");
% endif

## Step 1: Compute transformed and physical flux at solution points
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
    % for comp in range(ndims):
    % for phys in range(ndims):
        tF_upts[${upt}][${comp}][${var}] += smats_upts[${upt}][${nidx(comp,phys)}]*f[${phys}][${var}];
      % endfor
    % endfor
    ## % for phys in range(ndims):
    ##   f_upts[${upt}][${phys}][${var}] = f[${phys}][${var}];
    ## % endfor
  % endfor
}
% endfor

## Check
% if check:
% for upt in range(nupts):
  ## % for comp in range(ndims):
  ##   printf("smats ${upt} ${comp}_x=%.14e ${comp}_y=%.14e \n", smats_upts[${upt}][${nidx(comp,0)}], smats_upts[${upt}][${nidx(comp,1)}]);
  ## % endfor
  % for var in range(nvars):
    printf("tF_upts_${var} ${upt} = ${pnd} \n", ${','.join([f'tF_upts[{upt}][{comp}][{var}]' for comp in range(ndims)])});
  % endfor
% endfor
% endif

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

## Storage for N^D (discontinuous), N* (prescribed), and C^D (source in characteristic form) for all flux points
fpdtype_t N_D[${nfacefpts}][${nvars}];
fpdtype_t N_star[${nfacefpts}][${nvars}];
fpdtype_t C_D[${nfacefpts}][${nvars}];

## First pass: Compute N^D, N*, and C^D for all flux points
% for f, fpt_idx in enumerate(facefpts):
{
  int fidx = ${f};
<%
nE = bnorms[f,0]
nN = bnorms[f,1]
t1E = t1s[f,0]
t1N = t1s[f,1]
if ndims == 3:
    nC = bnorms[f,2]
    t1C = t1s[f,2]
    t2E = t2s[f,0]
    t2N = t2s[f,1]
    t2C = t2s[f,2]
%>

% if check:
  printf("\n Flux point %d\n", ${fpt_idx});
% endif

  ## Get face normals at flux point
  fpdtype_t norm_nl[${ndims}] = {${", ".join([f'normnl_ffpt[{f}][{i}]' for i in range(ndims)])}};
  fpdtype_t jac = jacs_ffpt[${f}];

  ## Check
% if check:
  printf("bnorm ${pnd} \n", ${','.join([f'{bnorms[f,comp]}' for comp in range(ndims)])});
  printf("t1 ${pnd} \n", ${','.join([f'{t1s[f,comp]}' for comp in range(ndims)])});
  % if ndims == 3:
  printf("t2 ${pnd} \n", ${','.join([f'{t2s[f,comp]}' for comp in range(ndims)])});
  % endif
  printf("norm_nl ${pnd} \n", ${','.join([f'norm_nl[{comp}]' for comp in range(ndims)])});
  printf("jac %.14e \n", jac);
% endif

  ## Check
  ## % if check:
  ##   % for upt in range(nupts):
    ## % for comp in range(ndims):
    ##   printf("smats_T ${upt} %.14e %.14e \n", smats_upts_T[${upt}][${nidx(comp,0)}], smats_upts_T[${upt}][${nidx(comp,1)}]);
    ## % endfor
  ##   % for var in range(nvars):
  ##     printf("tF_T ${upt} ${var} %.14e %.14e \n", tF_T[${upt}][0][${var}], tF_T[${upt}][1][${var}]);
  ##   % endfor
  ##   % endfor
  ## % endif


  ## Step 1a: Compute physical flux at our flux point
  fpdtype_t ul[${nvars}];
  fpdtype_t gradul[${ndims}][${nvars}];
  % for var in range(nvars):
    ul[${var}] = u_fpts[${fpt_idx}][${var}];
    ## TODO: Add later
    ## % for dim in range(ndims):
    ##   gradul[${dim}][${var}] = gradu_fpts[${nfidx(dim,nfpts)}][${var}];
    ## % endfor
  % endfor
  fpdtype_t fl[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'p', 'v')};
  ## TODO: Add later
  ## ${pyfr.expand('viscous_flux_add', 'ul', 'gradul', 'fl')};

  ## Check
% if check:
  ## % for var in range(nvars):
  ##   printf("grad_fpt ${fpt_idx} ${var} = %.2f %.2f \n",gradu_fpts[${nfidx(0,fpt_idx)}][${var}], gradu_fpts[${nfidx(1,fpt_idx)}][${var}]);
  ## % endfor

  ## Check
  % for var in range(nvars):
    printf("ul ${var} = %.14e \n", ul[${var}]);
    printf("fl ${var} = ${pnd} \n", ${','.join([f'fl[{phys}][{var}]' for phys in range(ndims)])});
  % endfor
    printf("p = %.14e  v = ${pnd} \n", p, ${','.join([f'v[{phys}]' for phys in range(ndims)])});
% endif

  ## Compute normal transformed flux
  fpdtype_t tfl_n[${nvars}] = {0};
  % for upt in range(nupts):
  % for var in range(nvars):
  % for comp in range(ndims):
    % if abs(m2[f,comp,upt]) > 0.0:
    tfl_n[${var}] += tF_upts[${upt}][${comp}][${var}]*${m2[f,comp,upt]};
    % endif
  % endfor
  % endfor
  % endfor

  ## Check
% if check:
  % for var in range(nvars):
    printf("tfl_n ${var} = %.14e \n", tfl_n[${var}]);
  % endfor
% endif

  ## Step 3: Compute derivative of face normal CS transformed flux, at flux point
  ## Also compute gradient of smats for geometric source term
  fpdtype_t dtFdE_full[${ndims}][${ndims}][${nvars}] = {{{0}}};
  fpdtype_t dsmatsdE_full[${ndims}][${ndims}][${ndims}] = {{{0}}};
  % for upt in range(nupts):
  {
    % for var in range(nvars):
    % for comp in range(ndims):
    % for comp2 in range(ndims):
      % if abs(m12[f,comp2,upt]) > 0.0:
        dtFdE_full[${comp}][${comp2}][${var}] += tF_upts[${upt}][${comp}][${var}]*${m12[f,comp2,upt]};
      % endif
    % endfor
    % endfor
    % endfor
    % for phys in range(ndims):
    % for comp in range(ndims):
    % for comp2 in range(ndims):
      % if abs(m12[f,comp2,upt]) > 0.0:
        dsmatsdE_full[${comp}][${comp2}][${phys}] += smats_upts[${upt}][${nidx(comp,phys)}]*${m12[f,comp2,upt]};
      % endif
    % endfor
    % endfor
    % endfor
  }
  % endfor

  ## Compute directional derivative to get derivative of normal transformed flux w.r.t. normal CS
  fpdtype_t dtFdE_T[${ndims}][${nvars}];
  % for var in range(nvars):
  {
    % if ndims == 2:
    dtFdE_T[0][${var}] = ${nE*nE}*   dtFdE_full[0][0][${var}] +
                         ${nN*nN}*   dtFdE_full[1][1][${var}] +
                         ${nE*nN}*  (dtFdE_full[0][1][${var}] +
                                     dtFdE_full[1][0][${var}]);
    dtFdE_T[1][${var}] = ${t1E*t1E}* dtFdE_full[0][0][${var}] +
                         ${t1N*t1N}* dtFdE_full[1][1][${var}] +
                         ${t1E*t1N}*(dtFdE_full[0][1][${var}] +
                                     dtFdE_full[1][0][${var}]);
    % else:
    dtFdE_T[0][${var}] = ${nE*nE}*   dtFdE_full[0][0][${var}] +
                         ${nN*nN}*   dtFdE_full[1][1][${var}] +
                         ${nC*nC}*   dtFdE_full[2][2][${var}] +
                         ${nE*nN}*  (dtFdE_full[0][1][${var}] + dtFdE_full[1][0][${var}]) +
                         ${nE*nC}*  (dtFdE_full[2][0][${var}] + dtFdE_full[0][2][${var}]) +
                         ${nN*nC}*  (dtFdE_full[2][1][${var}] + dtFdE_full[1][2][${var}]);
    dtFdE_T[1][${var}] = ${t1E*t1E}* dtFdE_full[0][0][${var}] +
                         ${t1N*t1N}* dtFdE_full[1][1][${var}] +
                         ${t1C*t1C}* dtFdE_full[2][2][${var}] +
                         ${t1E*t1N}*(dtFdE_full[0][1][${var}] + dtFdE_full[1][0][${var}]) +
                         ${t1E*t1C}*(dtFdE_full[2][0][${var}] + dtFdE_full[0][2][${var}]) +
                         ${t1N*t1C}*(dtFdE_full[2][1][${var}] + dtFdE_full[1][2][${var}]);
    dtFdE_T[2][${var}] = ${t2E*t2E}* dtFdE_full[0][0][${var}] +
                         ${t2N*t2N}* dtFdE_full[1][1][${var}] +
                         ${t2C*t2C}* dtFdE_full[2][2][${var}] +
                         ${t2E*t2N}*(dtFdE_full[0][1][${var}] + dtFdE_full[1][0][${var}]) +
                         ${t2E*t2C}*(dtFdE_full[2][0][${var}] + dtFdE_full[0][2][${var}]) +
                         ${t2N*t2C}*(dtFdE_full[2][1][${var}] + dtFdE_full[1][2][${var}]);

    % endif
  }
  % endfor

  fpdtype_t dsmatsdE_T[${ndims}];
  {
    % for phys in range(ndims):
      % if ndims == 2:
        dsmatsdE_T[${phys}] = ${nE*nE}* dsmatsdE_full[0][0][${phys}] +
                              ${nN*nN}* dsmatsdE_full[1][1][${phys}] +
                              ${nE*nN}*(dsmatsdE_full[0][1][${phys}] + dsmatsdE_full[1][0][${phys}]);
      % else:
        dsmatsdE_T[${phys}] = ${nE*nE}* dsmatsdE_full[0][0][${phys}] +
                              ${nN*nN}* dsmatsdE_full[1][1][${phys}] +
                              ${nC*nC}* dsmatsdE_full[2][2][${phys}] +
                              ${nE*nN}*(dsmatsdE_full[0][1][${phys}] + dsmatsdE_full[1][0][${phys}]) +
                              ${nE*nC}*(dsmatsdE_full[2][0][${phys}] + dsmatsdE_full[0][2][${phys}]) +
                              ${nN*nC}*(dsmatsdE_full[2][1][${phys}] + dsmatsdE_full[1][2][${phys}]);
      % endif
    % endfor
  }

  ## Check
% if check:
    % for var in range(nvars):
        printf("dtFdE_full ${var} 0 = ${pnd} \n", ${','.join([f'dtFdE_full[0][{comp}][{var}]' for comp in range(ndims)])});
        printf("dtFdE_full ${var} 1 = ${pnd} \n", ${','.join([f'dtFdE_full[1][{comp}][{var}]' for comp in range(ndims)])});
        % if ndims == 3:
        printf("dtFdE_full ${var} 2 = ${pnd} \n", ${','.join([f'dtFdE_full[2][{comp}][{var}]' for comp in range(ndims)])});
        % endif
    % endfor
  % for var in range(nvars):
    printf("dtFdE_T ${var} = ${pnd} \n", ${','.join([f'dtFdE_T[{comp}][{var}]' for comp in range(ndims)])});
  % endfor
  % for phys in range(ndims):
    printf("dsmatsdE_T ${phys} = %.14e \n", dsmatsdE_T[${phys}]);
  % endfor
% endif

  ## Step 4: Compute the characteristic wave strengths, N
  fpdtype_t N[${nvars}];
  fpdtype_t S[${nvars}];
  fpdtype_t source[${nvars}];

  fpdtype_t dEdE[${nvars}] = {${','.join([f'dtFdE_T[0][{i}]' for i in range(nvars)])}};
  fpdtype_t dGdN[${nvars}] = {${','.join(['+'.join([f'dtFdE_T[{j+1}][{i}]' for j in range(ndims-1)]) for i in range(nvars)])}};

  % for var in range(nvars):
  {
    source[${var}] = ${'+'.join([f'fl[{dim}][{var}]*dsmatsdE_T[{dim}]' for dim in range(ndims)])};
    dEdE[${var}] -= source[${var}];
    dGdN[${var}] += source[${var}];
  }
  % endfor

  ${pyfr.expand('WU_dot_dE','dEdE','N','ul','p','v')}
  ${pyfr.expand('WU_dot_dE','dGdN','S','ul','p','v')}

  ## Store N^D (before BC application) 
  % for var in range(nvars):
    N_D[${f}][${var}] = N[${var}];
  % endfor

  ## Compute and store C^D = WU·(tfl_n)
  fpdtype_t C_D_temp[${nvars}];
  ${pyfr.expand('WU_dot_dE','tfl_n','C_D_temp','ul','p','v')}
  % for var in range(nvars):
    C_D[${f}][${var}] = C_D_temp[${var}];
  % endfor

  ## Check
% if check:
  % for var in range(nvars):
    printf("source ${var} %.14e \n", source[${var}]);
  % endfor
  % for var in range(nvars):
    printf("N_D[${f}][${var}] = %.14e \n", N_D[${f}][${var}]);
    printf("C_D[${f}][${var}] = %.14e \n", C_D[${f}][${var}]);
  % endfor
% endif

  ## Step 5: Compute wave amplitudes for unknown waves (applies BC to get N*)
  ${pyfr.expand('compute_wave_amp', 'ul', 'p', 'v', 'jac', 'N', 'S', 'norm_nl')};

  ## Store N* (after BC application)
  % for var in range(nvars):
    N_star[${f}][${var}] = N[${var}];
  % endfor

  ## Check
% if check:
  % for var in range(nvars):
    printf("N_star[${f}][${var}] = %.14e \n", N_star[${f}][${var}]);
  % endfor
% endif

}
% endfor

## Second pass: Solve the coupled system for common normal flux
## F_fp^perp = WU^{-1} * [C^D + G^{-1} * (N* - N^D)]_fp

## Step 1: Compute G^{-1} * (N* - N^D) for each variable
fpdtype_t G_inv_dN[${nfacefpts}][${nvars}];
% for var in range(nvars):
  % for i in range(nfacefpts):
    G_inv_dN[${i}][${var}] = 0.0;
    % for j in range(nfacefpts):
      G_inv_dN[${i}][${var}] += ${G_inv[i,j]} * (N_star[${j}][${var}] - N_D[${j}][${var}]);
    % endfor
  % endfor
% endfor

## Step 2: For each flux point, compute the common normal flux
% for f, fpt_idx in enumerate(facefpts):
{
  ## Get face normals at flux point (needed for WU^{-1} transformation)
  fpdtype_t norm_nl[${ndims}] = {${", ".join([f'normnl_ffpt[{f}][{i}]' for i in range(ndims)])}};
  
  ## Get solution at flux point for transformation
  fpdtype_t ul[${nvars}];
  % for var in range(nvars):
    ul[${var}] = u_fpts[${fpt_idx}][${var}];
  % endfor
  fpdtype_t fl[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'p', 'v')};

  ## Compute characteristic amplitudes for this flux point
  fpdtype_t N_common[${nvars}];
  % for var in range(nvars):
    N_common[${var}] = C_D[${f}][${var}] + G_inv_dN[${f}][${var}];
  % endfor

  ## Transform back to conservative variables: F_common = WU^{-1} * N_common
  fpdtype_t F_common[${nvars}];
  ${pyfr.expand('WUinv_dot_N','N_common','F_common','ul','p','v')};

  ## Store the common normal flux in u_fpts (this is the output)
  % for var in range(nvars):
    u_fpts[${fpt_idx}][${var}] = F_common[${var}];
    % if check:
      printf("F_common[${fpt_idx}][${var}] = %.14e \n", F_common[${var}]);
    % endif
  % endfor
}
% endfor

</%pyfr:kernel>
