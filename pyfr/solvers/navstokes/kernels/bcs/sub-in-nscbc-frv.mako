<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>
<% sq2 = 2**0.5 %>
<% invsq2 = 2**-0.5 %>

## <% check = True %>

<%pyfr:macro name='compute_wave_amp' params='u, p, v, jac, N, S, norm_nl'>
  fpdtype_t nx = norm_nl[0];
  fpdtype_t ny = norm_nl[1];

  fpdtype_t rhoR = jac*${c['K_rho']}*(${c['rho']} - rho);
  fpdtype_t uR = jac*${c['K_u']}*(${c['u']} - v[0]);
  fpdtype_t vR = jac*${c['K_v']}*(${c['v']} - v[1]);

  ## Need to solve for all but outgoing wave
% if ndims == 2:
  % if check:
  printf("rhoR %e uR %e  vR %e\n", rhoR, uR, vR);
  % endif

  N[0] = (-(rhoR) + rho*invc*((uR*nx + vR*ny) - ${sq2}*(N[3])));
  N[1] = (uR*ny - vR*nx);
  N[2] = (N[3] - ${sq2}*(uR*nx + vR*ny));

% elif ndims == 3:

  fpdtype_t nz = norm_nl[2];
  fpdtype_t wR = jac*${c['K_w']}*(${c['w']} - v[2]);

  % if check:
  printf("rhoR %e uR %e  vR %e wR %e\n", rhoR, uR, vR, wR);
  % endif

  N[0] = -invc*(${sq2}*N[4]*nx*rho + (ny*ny+nz*nz-1.0)*rho*uR + c*(nx*rhoR + nz*vR - ny*wR) - nx*rho*(ny*vR + nz*wR));
  N[1] =  invc*(${-sq2}*N[4]*ny*rho - c*(ny*rhoR - nz*uR + nx*wR) + ny*rho*(nx*uR + ny*vR + nz*wR));
  N[2] =  invc*(${-sq2}*N[4]*nz*rho - c*(nz*rhoR + ny*uR - nx*vR) + nz*rho*(nx*uR + ny*vR + nz*wR));
  N[3] = N[4] - ${sq2}*(nx*uR + ny*vR + nz*wR);

% endif
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
  % for i in range(nvars):
    ur[${i}] = ul[${i}];
  % endfor
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>