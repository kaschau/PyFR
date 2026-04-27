<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.nscbc_fr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.nscbc_decomp'/>

<%def name="nscbc_body()">
## Compute transformed flux at solution points
fpdtype_t tf_upts[${nupts}][${ndims}][${nvars}] = {{{0}}};
% for upt in range(nupts):
{
  fpdtype_t u[${nvars}];
  % for var in range(nvars):
  u[${var}] = u_upts[${upt}][${var}];
  % endfor
  fpdtype_t gradu[${ndims}][${nvars}];
  % for var in range(nvars):
    % for dim in range(ndims):
  gradu[${dim}][${var}] = gradu_upts[${dim*nupts + upt}][${var}];
    % endfor
  % endfor
  fpdtype_t q[${nvars + 2}];
  fpdtype_t qh[${4 + mcf.ns}];
  ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};
  fpdtype_t f[${ndims}][${nvars}];
  ${pyfr.expand('inviscid_flux', 'u', 'f', 'q')};
  fpdtype_t qt[${mcf.ns + 2}];
  ${pyfr.expand('mixture_transport', 'u', 'q', 'qh', 'qt')};
  ${pyfr.expand('viscous_flux_add', 'u', 'gradu', 'q', 'qh', 'qt', 'f')};
  % for var in range(nvars):
    % for comp in range(ndims):
      % for phys in range(ndims):
  tf_upts[${upt}][${comp}][${var}] += smats_upts[${upt}][${comp*ndims + phys}]*f[${phys}][${var}];
      % endfor
    % endfor
  % endfor
}
% endfor

fpdtype_t tnf_D[${nfpts}][${nvars}] = {{0}};
${pyfr.expand('disc_nflux', 'tf_upts', 'tnf_D')}

fpdtype_t delta_div[${nfacefpts}][${nvars}];
% for f, fpt_idx in enumerate(facefpts):
{
  fpdtype_t jac = jacs_ffpts[${f}];

  fpdtype_t norm_nl[${ndims}];
  fpdtype_t t1[${ndims}];
% if ndims == 3:
  fpdtype_t t2[${ndims}];
% endif
  ${pyfr.expand('face_cs', 'smats_upts', 'norm_nl', 't1', 't2', f)}

  ## Load state at flux point and compute primitives via MC stateFrom-cons
  fpdtype_t ul[${nvars}];
  % for var in range(nvars):
  ul[${var}] = u_fpts[${fpt_idx}][${var}];
  % endfor
  fpdtype_t ql[${nvars + 2}];
  fpdtype_t qhl[${4 + mcf.ns}];
  ${pyfr.expand('stateFrom-cons', 'ul', 'ql', 'qhl')};

  fpdtype_t div[${nvars}];
  ${pyfr.expand('ref_div', 'tf_upts', 'div', f)}

  ## Project divergence onto characteristics, apply BC, project back
  fpdtype_t Phi[${nvars + 1}];
  ${pyfr.expand(f'WU_dot_div-{decomp_type}','div','Phi','ql','qhl')}

  ${pyfr.expand('compute_wave_amp', 'ul', 'ql', 'qhl', 'Phi', 'jac')};

  fpdtype_t div_star[${nvars}];
  ${pyfr.expand(f'WUinv_dot_Phi-{decomp_type}','Phi','div_star','ql','qhl')};

  % for var in range(nvars):
  delta_div[${f}][${var}] = div_star[${var}] - div[${var}];
  % endfor
}
% endfor

${pyfr.expand('fr_update', 'u_fpts', 'tnf_D', 'delta_div')}

</%def>
