<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>
<%from math import sqrt%>

<%pyfr:macro name='set_normal' params='bnorm, norm_nl'>
 ## For inlets, we want inward facing normals
 % for dim in range(ndims):
   bnorm[${dim}] *= -1.0;
   norm_nl[${dim}] *= -1.0;
 % endfor
</%pyfr:macro>

<%pyfr:macro name='compute_wave_amp' params='u, p, v, jac, N, S, norm_nl'>
  fpdtype_t nx = norm_nl[0];
  fpdtype_t ny = norm_nl[1];

  fpdtype_t rhoR = jac*${c['sigma']}*(${c['rho']} - rho);
  fpdtype_t uR = jac*${c['sigma']}*(${c['u']} - v[0]);
  fpdtype_t vR = jac*${c['sigma']}*(${c['v']} - v[1]);

  ## Need to solve for all but outgoing wave
  N[0] = -(rhoR+S[0]) - rho*invc*(${sqrt(2)}*(N[3]+S[3]) - (uR*nx + vR*ny));
  N[1] = -vR*nx + uR*ny - S[1];
  N[2] = N[3] + S[3] - ${sqrt(2)}*(uR*nx+vR*ny) - S[2];

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
  % for i in range(nvars):
    ur[${i}] = ul[${i}];
  % endfor
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
