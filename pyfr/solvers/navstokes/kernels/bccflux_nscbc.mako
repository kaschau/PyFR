<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>
<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<%include file='pyfr.solvers.navstokes.kernels.bcs.${bctype}'/>

<%def name="nidx(comp,phys)">
  <% return comp*ndims + phys %>
</%def>\
<%def name="nvidx(dim,fpt)">
  <% return dim*nfpts + fpt %>
</%def>\
<% from math import sqrt %>\


<%pyfr:macro name='mmmul' params='A, B, C'>
% for row in range(nvars):
% for col in range(nvars):
  C[${row}][${col}] = ${'+'.join([f'A[{row}][{c}]*B[{c}][{col}]' for c in range(nvars)])};
% endfor
% endfor
</%pyfr:macro>

<%pyfr:macro name='mvmul' params='A, b, c'>
% for row in range(nvars):
% for col in range(nvars):
c[${row}] += A[${row}][${col}]*b[${col}];
% endfor
% endfor
</%pyfr:macro>

<%pyfr:kernel name='bccflux_nscbc' ndim='1'
              u_ele='in view fpdtype_t[${str(nupts)}][${str(nvars)}]'
              u_fpt='inout view fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              gradu_fpt='in view fpdtype_t[${str(ndims*nfpts)}][${str(nvars)}]'
              nl_ffpt='in fpdtype_t[${str(nfacefpts)}][${str(ndims)}]'
              smats_ele='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              jacs_ffpt='in fpdtype_t[${str(nfacefpts)}]'>

## printf("\n*************ELEMENT************\n");

## Step 1: Compute transformed visc flux, pressure, velocity at solution points
fpdtype_t tFi[${nupts}][${ndims}][${nvars}] = {{{0}}};
for (int uidx = 0; uidx < ${nupts}; uidx++)
{
  fpdtype_t ul[${nvars}];
  % for var in range(nvars):
    ul[${var}] = u_ele[uidx][${var}];
  % endfor
  fpdtype_t fi[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  ${pyfr.expand('inviscid_flux', 'ul', 'fi', 'p', 'v')};

  % for var in range(nvars):
    % for comp in range(ndims):
      % for phys in range(ndims):
        tFi[uidx][${comp}][${var}] += smats_ele[uidx][${nidx(comp,phys)}]*fi[${phys}][${var}];
      % endfor
    % endfor
  % endfor
}

## Check
## % for upt in range(nupts):
##   % for var in range(nvars):
##     printf("upt ${upt} tFi_${var} = %.14e %.14e \n", tFi[${upt}][0][${var}], tFi[${upt}][1][${var}]);
##   % endfor
## % endfor

## % for upt in range(nupts):
## % for var in range(nvars):
##   printf("upt ${upt} ${var} = %.2f \n", u_ele[${upt}][${var}]);
## % endfor
## % endfor

## Iterate over the flux points on our face
% for f, fpt_idx in enumerate(facefpts):
{

  ## printf("\n Flux point %d\n", ${fpt_idx});

  ## Step 1a: Compute physical flux point values
  fpdtype_t fl[${ndims}][${nvars}];
  fpdtype_t p, v[${ndims}];
  fpdtype_t ul[${nvars}];
  % for var in range(nvars):
    ul[${var}] = u_fpt[${fpt_idx}][${var}];
  % endfor
  ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'p', 'v')};

  ## Check
  ## % for var in range(nvars):
  ##   printf("ul ${var} = %f\n", ul[${var}]);
  ##   printf("fl ${var} = %f %f \n", fl[0][${var}], fl[1][${var}]);
  ## % endfor
  ##   printf("p = %f  v = %f %f \n", p, v[0], v[1]);

  fpdtype_t nl[${ndims}] = {${", ".join([f'nl_ffpt[{f}][{i}]' for i in range(ndims)])}};
  fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
  fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

  fpdtype_t c = sqrt(${c['gamma']}*p/ul[0]);
  fpdtype_t rho = ul[0];
  fpdtype_t invrho = 1.0/rho;
  fpdtype_t invc = 1.0/c;
  fpdtype_t invcsq = invc*invc;

  fpdtype_t dQdU[${nvars}][${nvars}] = {
    {1.0, 0.0, 0.0, 0.0},
    {-v[0]*invrho, invrho, 0.0, 0.0},
    {-v[1]*invrho, 0.0, invrho, 0.0},
    {${(c['gamma']-1.0)/2.0}*(${pyfr.dot('v[{i}]', i=ndims)}), ${1.0-c['gamma']}*v[0], ${1.0-c['gamma']}*v[1], ${c['gamma']-1.0}},
  };
  fpdtype_t dUdQ[${nvars}][${nvars}] = {
    {1.0, 0.0, 0.0, 0.0},
    {v[0], rho, 0.0, 0.0},
    {v[1], 0.0, rho, 0.0},
    {0.5*(${pyfr.dot('v[{i}]', i=ndims)}), ul[0], ul[1], ${1.0/(c['gamma']-1)}}
  };
  fpdtype_t dWdQ[${nvars}][${nvars}] = {
    {norm_nl[0], 0.0, 0.0, -norm_nl[0]*invcsq},
    {0.0, -norm_nl[1], norm_nl[0], 0.0},
    {0.0,  norm_nl[0]*${1.0/sqrt(2)},  norm_nl[1]*${1.0/sqrt(2)}, invrho*invc*${1.0/sqrt(2)}},
    {0.0, -norm_nl[0]*${1.0/sqrt(2)}, -norm_nl[1]*${1.0/sqrt(2)}, invrho*invc*${1.0/sqrt(2)}},
  };
  fpdtype_t dQdW[${nvars}][${nvars}] = {
    {norm_nl[0], 0.0, rho*invc*${1.0/sqrt(2)}, rho*invc*${1.0/sqrt(2)}},
    {0.0, -norm_nl[1], norm_nl[0]*${1.0/sqrt(2)}, -norm_nl[0]*${1.0/sqrt(2)}},
    {0.0,  norm_nl[0], norm_nl[1]*${1.0/sqrt(2)}, -norm_nl[1]*${1.0/sqrt(2)}},
    {0.0, 0.0, rho*c*${1.0/sqrt(2)}, rho*c*${1.0/sqrt(2)}},
  };

  ## Check
  ## % for row in range(nvars):
  ##   printf("WQ[${row}] = %f %f %f %f \n", dWdQ[${row}][0],dWdQ[${row}][1],dWdQ[${row}][2],dWdQ[${row}][3]);
  ## % endfor
  ## % for row in range(nvars):
  ##   printf("QU[${row}] = %f %f %f %f \n", dQdU[${row}][0],dQdU[${row}][1],dQdU[${row}][2],dQdU[${row}][3]);
  ## % endfor

  fpdtype_t PU[${nvars}][${nvars}] = {{0}};
  ${pyfr.expand('mmmul','dWdQ','dQdU','PU')};
  fpdtype_t PUinv[${nvars}][${nvars}] = {{0}};
  ${pyfr.expand('mmmul','dUdQ','dQdW','PUinv')};

  ## Check
  ## % for row in range(nvars):
  ##   printf("PU[${row}] = %f %f %f %f \n", PU[${row}][0],PU[${row}][1],PU[${row}][2],PU[${row}][3]);
  ## % endfor

  fpdtype_t fl_n[${nvars}] = {0};
  % for upt in range(nupts):
  % for comp in range(ndims):
  % for var in range(nvars):
    fl_n[${var}] += tFi[${upt}][${comp}][${var}]*${m2[f,comp,upt]};
  % endfor
  % endfor
  % endfor

  ## Check
  ## % for var in range(nvars):
  ##   printf("fl ${var} = %f %f \n", fl[0][${var}], fl[1][${var}]);
  ##   printf("fl_n ${var} = %f \n", fl_n[${var}]);
  ## % endfor
  ##   printf("p = %f  v = %f %f \n", p, v[0], v[1]);

  ## ## Step 2: Compute derivative in transformed space of tflux, at flux point
  ## Also compute gradient of smats for geometric source term
  fpdtype_t dtFidE[${ndims}][${nvars}] = {{0}};
  fpdtype_t dsmatsdE[${ndims}][${ndims}] = {{0}};
  % for upt in range(nupts):
    % for comp in range(ndims):
      % for var in range(nvars):
        dtFidE[${comp}][${var}] += tFi[${upt}][${comp}][${var}]*${m12[f, comp, upt]};
      % endfor
      % for phys in range(ndims):
        dsmatsdE[${comp}][${phys}] += smats_ele[${upt}][${nidx(comp,phys)}]*${m12[f, comp, upt]};
      % endfor
    % endfor
  % endfor

  ## Check
  ## % for var in range(nvars):
  ##   printf("dtFidE_${var} = %f %f \n", dtFidE[0][${var}], dtFidE[1][${var}]);
  ## % endfor

  ## ## Step 3: Realign derivatives to a face-normal orientation where \Xi is normal to face
  fpdtype_t bnorm[${ndims}] = {${','.join([str(i) for i in bnorm_facefpts[f,:]])}};
  fpdtype_t dtFidE_n[${ndims}][${nvars}];
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

  ## Check
  ## % for var in range(nvars):
  ##   printf(" dtFidE_n ${var} %.14e %.14e\n", dtFidE_n[0][${var}], dtFidE_n[1][${var}]);
  ## % endfor

  ## fpdtype_t dsmatsdE_n[${ndims}][${nvars}];
  ## % for phys in range(ndims):
  ## {
  ##   fpdtype_t dsmats_temp[${ndims}] = {${','.join([f'dsmatsdE[{comp}][{phys}]' for comp in range(ndims)])}};
  ##   fpdtype_t dsmats_n_temp[${ndims}];
  ##   ${pyfr.expand('transform_to', 'bnorm', 'dsmats_temp', 'dsmats_n_temp', off=0)};
  ##   % for comp in range(ndims):
  ##     dsmatsdE_n[${comp}][${phys}] = dsmats_n_temp[${comp}];
  ##   % endfor
  ## }
  ## % endfor

  ## ## Step 4: Compute the characteristic wave strengths, N

  fpdtype_t N[${nvars}] = {0};
  fpdtype_t S[${nvars}] = {0};
  {
    fpdtype_t dEdE[${nvars}] = {${','.join([f'dtFidE_n[0][{i}]' for i in range(nvars)])}};
    ${pyfr.expand('mvmul','PU','dEdE','N')}
    fpdtype_t dGdN[${nvars}] = {${','.join([f'dtFidE_n[1][{i}]' for i in range(nvars)])}};
    % if ndims == 3:
      fpdtype_t dGdN[${nvars}] = {${','.join([f'dtFidE_n[2][{i}]' for i in range(nvars)])}};
    % endif
    ${pyfr.expand('mvmul','PU','dGdN','S')}
  }

  ## Step 5: Replace incoming wave amplitudes
  ## ${pyfr.expand('compute_L', 'L', 'u_f')};
  fpdtype_t Msq = (${pyfr.dot('v[{i}]', i=ndims)})*invcsq;
  N[${nvars-1}] = jacs_ffpt[${f}]*${c['sigma']/sqrt(2)}/ul[0]*(1.0-Msq)*(p - ${c['p']}) - S[${nvars-1}];

  ## Check
  ## % for i in range(nvars):
  ##   printf("N_${i} = %.14e \n", N[${i}]);
  ## % endfor

  ## ## Step 6: Compute dFstardE values normal to face
  {
    fpdtype_t dEdE[${nvars}] = {0};
    ${pyfr.expand('mvmul','PUinv','N','dEdE')};
    % for var in range(nvars):
      dtFidE_n[0][${var}] = dEdE[${var}];
    % endfor
  }


  ## ## Convert back
  fpdtype_t dtFidE_star[${ndims}][${nvars}];
  % for var in range(nvars):
  {
    ## Correct signs
    fpdtype_t dF_n_temp[${ndims}] = {${','.join([f'dtFidE_n[{i}][{var}]' for i in range(ndims)])}};
    fpdtype_t dF_temp[${ndims}];
    ${pyfr.expand('transform_from', 'bnorm', 'dF_n_temp', 'dF_temp', off=0)};
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
    u_fpt[${fpt_idx}][${var}] = ${'+'.join([f'dtFidE_star[{dim}][{var}]' for dim in range(ndims)])};
    ## printf("dtFidE_* %.14e \n", ${'+'.join([f'dtFidE_star[{dim}][{var}]' for dim in range(ndims)])});

    ## Add geomteric source term
    ## uf[${fpt_idx}][${var}] += (${pyfr.dot('f_f[{i}][var]','dsmatsdE_n[0][{i}]', i=ndims)});
    ## printf("geom source %f\n",${pyfr.dot('f_f[{i}][var]','dsmatsdE_n[0][{i}]', i=ndims)});

    ## subtract our flux gradient on the face (from interior values)
    u_fpt[${fpt_idx}][${var}] -= (${'+'.join([f'dtFidE[{dim}][{var}]' for dim in range(ndims)])});
    ## printf("dtFidE %f \n", (${'+'.join([f'dtFidE[{dim}][{var}]' for dim in range(ndims)])}));

    ## divide by ~\del \dot g
    u_fpt[${fpt_idx}][${var}] /= ${m11[f]};

    ## add the transformed, normal flux from interior values
    ## fpdtype_t fl_n2 = ${pyfr.dot('nl[{i}]', 'fl[{i}][var]', i=ndims)};
    ## u_fpt[${fpt_idx}][${var}] += fl_n2;
    u_fpt[${fpt_idx}][${var}] += fl_n[${var}];

    ## Check
    ## printf("f_n =  %.14e \n", u_fpt[${fpt_idx}][${var}]);
  }
  % endfor
}
% endfor

</%pyfr:kernel>
