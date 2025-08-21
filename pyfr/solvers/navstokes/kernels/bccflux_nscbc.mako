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
              normnl_ffpts='in fpdtype_t[${str(nfacefpts)}][${str(ndims)}]'
              smats_upts='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              jacs_ffpts='in fpdtype_t[${str(nfacefpts)}]'>

% if check:
printf("\n*************ELEMENT************\n");
% endif

## Step 1: Compute transformed flux at solution points
fpdtype_t tF_upts[${nupts}][${ndims}][${nvars}] = {{{0}}};
% for upt in range(nupts):
{
  fpdtype_t u[${nvars}];
  fpdtype_t gradu[${ndims}][${nvars}];
  % for var in range(nvars):
    u[${var}] = u_upts[${upt}][${var}];
    ## TODO
    ## % for dim in range(ndims):
    ##   gradu[${dim}][${var}] = gradu_upts[${nuidx(dim,upt)}][${var}];
    ## % endfor
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
  % endfor
}
% endfor

## Check
% if check:
% for upt in range(nupts):
  % for var in range(nvars):
    printf("tF_upts_${var} ${upt} = ${pnd} \n", ${','.join([f'tF_upts[{upt}][{comp}][{var}]' for comp in range(ndims)])});
  % endfor
% endfor
% endif

## Step 2: Compute transformed, discontinuous, normal flux at all flux points
fpdtype_t tf_TD[${nfpts}][${nvars}] = {{0}};
% for fpt_idx in fptidx:
{
  ## Compute normal transformed flux
  % for upt in range(nupts):
  % for var in range(nvars):
  % for comp in range(ndims):
    % if abs(m2[fpt_idx,comp,upt]) > 0.0:
    tf_TD[${fpt_idx}][${var}] += tF_upts[${upt}][${comp}][${var}]*${m2[fpt_idx,comp,upt]};
    % endif
  % endfor
  % endfor
  % endfor
}
% endfor

## Check
% if check:
% for f in range(nfpts):
  printf("\n Total Flux point %d\n", ${f});
  % for var in range(nvars):
    printf("tf_TD ${var} = %.14e \n", tf_TD[${f}][${var}]);
  % endfor
% endfor
% endif

## Step 3: Compute NSCBC boundary corrections at face flux points
fpdtype_t ddtF_TdE[${nfacefpts}][${nvars}];
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
  printf("\n Face Flux point %d\n", ${fpt_idx});
% endif

  ## Get face normals at flux point
  fpdtype_t norm_nl[${ndims}] = {${", ".join([f'normnl_ffpts[{f}][{i}]' for i in range(ndims)])}};
  fpdtype_t jac = jacs_ffpts[${f}];

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

  ## Compute physical flux at our flux point
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

  ## Compute derivative of face normal CS transformed flux, at flux point
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

  ## Step 5: Compute wave amplitudes for unknown waves (applies BC to get N*)
  ${pyfr.expand('compute_wave_amp', 'ul', 'p', 'v', 'jac', 'N', 'S', 'norm_nl')};

  ## Check
% if check:
  % for var in range(nvars):
    printf("N_star var=${var} = %.14e \n", N[${var}]);
  % endfor
% endif

  ## Step 6: Compute dtFdE_T* values normal to face
  fpdtype_t dtFdE_Ts[${nvars}];
  {
    ${pyfr.expand('WUinv_dot_N','N','dtFdE_Ts','ul','p','v')};
    % for var in range(nvars):
      dtFdE_Ts[${var}] += source[${var}];
    % endfor
  }

  ## Now solve for derivatives in the primary CS
  fpdtype_t dtFdE[${ndims}][${nvars}] = {{0}};
    % for var in range(nvars):
      % if ndims == 2:
        {
          fpdtype_t dFstar = dtFdE_Ts[${var}];
          fpdtype_t dFdN = dtFdE_full[0][1][${var}];
          fpdtype_t dGdE = dtFdE_full[1][0][${var}];
          fpdtype_t dGdN_T = dtFdE_T[1][${var}];
          dtFdE[0][${var}] = -((dFstar - (dFdN + dGdE)*${nE*nN})*${t1N*t1N}-${nN*nN}*(-((dFdN+dGdE)*${t1E*t1N})+dGdN_T))/${nN*nN*t1E*t1E-nE*nE*t1N*t1N};
          dtFdE[1][${var}] = (${-t1E}*(dFstar*${t1E} + (dFdN+dGdE)*${nE*(-nN*t1E+nE*t1N)})+${nE*nE}*dGdN_T)/${-nN*nN*t1E*t1E+nE*nE*t1N*t1N};
        }
      % else:
        {
          fpdtype_t dFstar = dtFdE_Ts[${var}];
          fpdtype_t dFdN = dtFdE_full[0][1][${var}];
          fpdtype_t dFdC = dtFdE_full[0][2][${var}];
          fpdtype_t dGdE = dtFdE_full[1][0][${var}];
          fpdtype_t dGdC = dtFdE_full[1][2][${var}];
          fpdtype_t dHdE = dtFdE_full[2][0][${var}];
          fpdtype_t dHdN = dtFdE_full[2][1][${var}];
          fpdtype_t dGdN_T = dtFdE_T[1][${var}];
          fpdtype_t dHdC_T = dtFdE_T[2][${var}];

          fpdtype_t facn  = ${nE*nN}  *(dGdE+dFdN) + ${nE*nC}  *(dHdE+dFdC) + ${nC*nN}  *(dGdC + dHdN);
          fpdtype_t fact1 = ${t1E*t1N}*(dGdE+dFdN) + ${t1E*t1C}*(dHdE+dFdC) + ${t1C*t1N}*(dGdC + dHdN);
          fpdtype_t fact2 = ${t2E*t2N}*(dGdE+dFdN) + ${t2E*t2C}*(dHdE+dFdC) + ${t2C*t2N}*(dGdC + dHdN);

          dtFdE[0][${var}] = (-fact2*${nN*nN*t1C*t1C} +
                               fact2*${nC*nC*t1N*t1N} +
                               fact1*${nN*nN*t2C*t2C} +
                               dFstar*${t1N*t1N*t2C*t2C} -
                               facn  *${t1N*t1N*t2C*t2C} -
                               fact1 *${nC*nC*t2N*t2N} -
                               dFstar*${t1C*t1C*t2N*t2N} +
                               facn  *${t1C*t1C*t2N*t2N} +
                              ${-nN*nN*t2C*t2C + nC*nC*t2N*t2N}*dGdN_T + ${nN*nN*t1C*t1C - nC*nC*t1N*t1N}*dHdC_T)/
                              ${nN*nN*(-t1E*t1E*t2C*t2C + t1C*t1C*t2E*t2E) + \
                                nE*nE*( t1N*t1N*t2C*t2C - t1C*t1C*t2N*t2N) + \
                                nC*nC*(-t1N*t1N*t2E*t2E + t1E*t1E*t2N*t2N)};

          dtFdE[1][${var}] = (fact2*${nE*nE*t1C*t1C} -
                              fact2*${nC*nC*t1E*t1E} -
                              fact1*${nE*nE*t2C*t2C} -
                              dFstar*${t1E*t1E*t2C*t2C} +
                              facn  *${t1E*t1E*t2C*t2C} +
                              fact1 *${nC*nC*t2E*t2E} +
                              dFstar*${t1C*t1C*t2E*t2E} -
                              facn  *${t1C*t1C*t2E*t2E} +
                             ${nE*nE*t2C*t2C - nC*nC*t2E*t2E}*dGdN_T + ${-nE*nE*t1C*t1C + nC*nC*t1E*t1E}*dHdC_T)/
                             ${nN*nN*(-t1E*t1E*t2C*t2C + t1C*t1C*t2E*t2E) + \
                               nE*nE*( t1N*t1N*t2C*t2C - t1C*t1C*t2N*t2N) + \
                               nC*nC*(-t1N*t1N*t2E*t2E + t1E*t1E*t2N*t2N)};

          dtFdE[2][${var}] = (fact2*${nN*nN*t1E*t1E} -
                              fact2*${nE*nE*t1N*t1N} -
                              fact1*${nN*nN*t2E*t2E} -
                              dFstar*${t1N*t1N*t2E*t2E} +
                              facn  *${t1N*t1N*t2E*t2E} +
                              fact1 *${nE*nE*t2N*t2N} +
                              dFstar*${t1E*t1E*t2N*t2N} -
                              facn  *${t1E*t1E*t2N*t2N} +
                             ${nN*nN*t2E*t2E - nE*nE*t2N*t2N}*dGdN_T + ${-nN*nN*t1E*t1E + nE*nE*t1N*t1N}*dHdC_T) /
                             ${nN*nN*(-t1E*t1E*t2C*t2C + t1C*t1C*t2E*t2E) + \
                               nE*nE*( t1N*t1N*t2C*t2C - t1C*t1C*t2N*t2N) + \
                               nC*nC*(-t1N*t1N*t2E*t2E + t1E*t1E*t2N*t2N)};
          }
      % endif
    % endfor

  % for var in range(nvars):
    ## we have dudt (~\del \dot ~f) at our flux point
    ddtF_TdE[${f}][${var}] = ${'+'.join([f'dtFdE[{comp}][{var}]' for comp in range(ndims)])};
    % if check:
        printf("dtFidE_* var=${var} %.14e \n", ddtF_TdE[${f}][${var}]);
    % endif
    ## subtract our flux gradient on the face (from discontinuous value)
    ddtF_TdE[${f}][${var}] -= ${'+'.join([f'dtFdE_full[{comp}][{comp}][{var}]' for comp in range(ndims)])};

    % if check:
        printf("dtFidE %.14e \n", ${'+'.join([f'dtFdE_full[{comp}][{comp}][{var}]' for comp in range(ndims)])});
    % endif
  % endfor
}
% endfor

## PASS 2
% if check:
    printf("\n PASS 2 \n");
% endif
% for f,fpt_idx in enumerate(facefpts):
{
  % if check:
    printf("\n Flux point %d\n", ${fpt_idx});
  % endif

  ## Compute R = GB_inv * GI * (F^\perp - F^D)
  % for var in range(nvars):
  {
    % if check:
        printf("\n VAR ${var} \n");
    % endif

    ## Compute R = GB_inv * GI * \Delta f_interior
    fpdtype_t R = 0;
    % for j, intfpt_idx in enumerate(intfpts):
      % if abs(GB_inv_GI[f,j]) > 0.0:
      R += ${GB_inv_GI[f,j]}*(u_fpts[${intfpt_idx}][${var}] - tf_TD[${intfpt_idx}][${var}]);
      % endif
    % endfor
    % if check:
      printf("R = %.14e \n", R);
    % endif
    ## Compute A = GB_inv * d \nabla f
    fpdtype_t A = 0.0;
    % for j in range(nfacefpts):
      % if abs(GB_inv[f,j]) > 0.0:
      A += ${GB_inv[f,j]} * ddtF_TdE[${j}][${var}];
      % endif
    % endfor
    % if check:
      printf("A = %.14e \n", A);
    % endif

    ## Store the common normal flux in u_fpts (this is the output)
    u_fpts[${fpt_idx}][${var}] = tf_TD[${fpt_idx}][${var}] + A - R;

    ## Check
    % if check:
        printf("final flux =  %.14e \n", u_fpts[${fpt_idx}][${var}]);
    % endif
  }
  % endfor
}
% endfor

</%pyfr:kernel>
