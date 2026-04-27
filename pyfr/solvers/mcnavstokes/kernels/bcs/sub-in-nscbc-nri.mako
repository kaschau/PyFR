<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.common'/>
<% sq2 = 2**0.5 %>
<% invsq2 = 2**-0.5 %>
<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='compute_wave_amp' params='u, q, qh, Phi, jac' externs='ploc, t'>
  fpdtype_t rho = q[${rhoix}];
  fpdtype_t T = q[${Tix}];
  fpdtype_t c = qh[2];
  fpdtype_t invT = 1.0/T;
  fpdtype_t gmo = qh[0] - 1.0;
  fpdtype_t v[${ndims}] = {${','.join([f'q[{vix+i}]' for i in range(ndims)])}};
  fpdtype_t un = ${pyfr.dot('norm_nl[{i}]','v[{i}]', i=ndims)};
  fpdtype_t ut1 = ${pyfr.dot('t1[{i}]','v[{i}]', i=ndims)};

  ## Targets with acoustic forcing
  fpdtype_t T_target = ${c['T']} + T*gmo/c*${c['u_a']};
  fpdtype_t un_target = ${c['un']} + ${c['u_a']} + ${c['u_v']};

  fpdtype_t unR = ${c['K_ac']}*(un_target - un);
  unR += (2.0*${c['du_a_dt']} + ${c['du_v_dt']});
  fpdtype_t TR = ${c['K_ac']}*(T_target - T);

  ## Isentropic inlet: Phi_0 = 0
  Phi[0] = 0.0;

  ## Tangential relaxation
  fpdtype_t ut1R = jac*${c['K_ut']}*(0.0 - ut1);
  Phi[1] = -ut1R;
% if ndims == 3:
  fpdtype_t ut2 = ${pyfr.dot('t2[{i}]','v[{i}]', i=ndims)};
  fpdtype_t ut2R = jac*${c['K_ut']}*(0.0 - ut2);
  Phi[2] = -ut2R;
% endif

  ## Species waves and accumulate weighted sum for isentropic correction
  fpdtype_t sumYR = 0.0;
% for n, spn in enumerate(mcf.sp_names):
  {
    fpdtype_t YR = ${c['K_ac']}*(${c[spn]} - q[${n}]);
    Phi[${ndims + n}] = -jac*YR;
    sumYR += ${1.0/mcf[n].MW}*YR;
  }
% endfor

  ## Compute mixture molecular weight
  fpdtype_t inv_MWmix = ${'+'.join([f'q[{n}]*{1.0/mcf[n].MW}' for n in range(mcf.ns)])};
  fpdtype_t MWmix = 1.0/inv_MWmix;

  ## Combined isentropic form
  Phi[${nvars - 1}] = ${-invsq2}*jac*(unR + c/gmo*(invT*TR + MWmix*sumYR));
</%pyfr:macro>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>
  % for i in range(nvars):
    ur[${i}] = ul[${i}];
  % endfor
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
