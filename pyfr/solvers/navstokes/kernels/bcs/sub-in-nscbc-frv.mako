<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>
<%from math import sqrt%>

<%pyfr:macro name='compute_wave_amp' params='u, p, v, jac, JN, JS, norm_nl'>
  fpdtype_t nx = norm_nl[0];
  fpdtype_t ny = norm_nl[1];

  fpdtype_t rhoR = jac*${c['K_rho']}*(${c['rho']} - rho);
  fpdtype_t uR = jac*${c['K_u']}*(${c['u']} - v[0]);
  fpdtype_t vR = jac*${c['K_v']}*(${c['v']} - v[1]);

  ## Need to solve for all but outgoing wave
  JN[0] = (-(rhoR+jac*JS[0]) + rho*invc*((uR*nx + vR*ny) - jac*${sqrt(2)}*(JN[3] + JS[3])));
  JN[1] = (uR*ny - vR*nx - jac*JS[1]);
  JN[2] = jac*(JN[3] + JS[3]) - ${sqrt(2)}*(uR*nx + vR*ny) - jac*JS[2];
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
  % for i in range(nvars):
    ur[${i}] = ul[${i}];
  % endfor
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
