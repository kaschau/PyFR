<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>
<% sq2 = 2**0.5 %>
<% invsq2 = 2**-0.5 %>
<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='compute_wave_amp' params='u, q, qh, Phi, jac' externs='ploc, t'>
  fpdtype_t invc = 1.0/qh[2];
  fpdtype_t invrho = 1.0/q[${rhoix}];

  fpdtype_t pR = jac*${c['K_p']}*(${c['p']} - q[${pix}]);

  ## Only need incoming wave
  Phi[${nvars}] = ${-sq2}*invc*invrho*pR;

</%pyfr:macro>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
  % for i in range(nvars):
    ur[${i}] = ul[${i}];
  % endfor
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
