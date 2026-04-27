<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>
<% sq2 = 2**0.5 %>
<% invsq2 = 2**-0.5 %>
<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='compute_wave_amp' params='u, q, qh, Phi, jac' externs='ploc, t'>
  fpdtype_t nx = norm_nl[0];
  fpdtype_t ny = norm_nl[1];

  fpdtype_t rho = q[${rhoix}];
  fpdtype_t T = q[${Tix}];
  fpdtype_t invT = 1.0/T;
  fpdtype_t c = qh[2];
  fpdtype_t invc = 1.0/c;
  fpdtype_t gmo = qh[0] - 1.0;

  fpdtype_t TR = jac*${c['K_T']}*(${c['T']} - T);
  fpdtype_t uR = jac*${c['K_u']}*(${c['u']} - q[${vix+0}]);
  fpdtype_t vR = jac*${c['K_v']}*(${c['v']} - q[${vix+1}]);

  fpdtype_t invRmix = 1.0/(qh[1] - qh[1]/qh[0]);
  fpdtype_t sumPhiYiRi = 0.0;
  % for n, spn in enumerate(mcf.sp_names):
  {
    fpdtype_t YR = ${c['K_Y']}*(${c[spn]} - q[${n}]);
    sumPhiYiRi += jac*YR*${mcf.Ru/mcf[n].MW};
    Phi[${ndims+n}] = -jac*YR;
  }
  % endfor

  ## Need to solve for all but outgoing wave
% if ndims == 2:

  fpdtype_t veloR = nx*uR + ny*vR;
  Phi[0] = rho*(TR*invT + sumPhiYiRi*invRmix - invc*gmo*veloR);
  Phi[1] = uR*ny - vR*nx;
  Phi[${nvars-1}] = ${-sq2}*veloR;

% elif ndims == 3:

  fpdtype_t nz = norm_nl[2];
  fpdtype_t wR = jac*${c['K_w']}*(${c['w']} - q[${vix+2}]);

  fpdtype_t veloR = nx*uR + ny*vR + nz*wR;
  fpdtype_t temp = sumPhiYiRi*invRmix + TR*invT - invc*gmo*veloR;
  Phi[0] = nx*rho*(temp) - (nz*vR - ny*wR);
  Phi[1] = ny*rho*(temp) - (-nz*uR + nx*wR);
  Phi[2] = nz*rho*(temp) - (ny*uR - nx*vR);
  Phi[${nvars-1}] = ${-sq2}*veloR;

% endif
</%pyfr:macro>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
  % for i in range(nvars):
    ur[${i}] = ul[${i}];
  % endfor
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
