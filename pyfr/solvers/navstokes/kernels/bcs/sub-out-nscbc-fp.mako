<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>
<%from math import sqrt%>

<%pyfr:macro name='compute_wave_amp' params='u, p, v, jac, N, S, norm_nl'>
  fpdtype_t rho = u[0];
  fpdtype_t csq = ${c['gamma']}*p/rho;
  fpdtype_t Msq = (${pyfr.dot('v[{i}]', i=ndims)})*invcsq;
  fpdtype_t alpha = sqrt(Msq);

  fpdtype_t pR = jac*${c['K_p']/sqrt(2)}/u[0]*(1.0-Msq)*(p - ${c['p']});

  ## Only need incoming wave
  N[${nvars-1}] = pR - (1.0 - alpha)*jac*S[${nvars-1}];

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
  % for i in range(nvars):
    ur[${i}] = ul[${i}];
  % endfor
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
