<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%import math %>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>
<%
    ns, nvars = mcf.ns, mcf.ns + ndims + 1
    Rsum = ' + '.join(f'q[{n}]*{mcf.Ru/mcf[n].MW}' for n in range(ns))
    kesum = ' + '.join(f'q[{vix + i}]*q[{vix + i}]' for i in range(ndims))
    mws = [sp.MW for sp in mcf.species]
    mwinvs = [1.0/sp.MW for sp in mcf.species]
%>

## Analytic Jacobian of the chemical source term, jac[r*nvars + c] =
## d(MW_r*omega_r)/du_c in conserved variables u = (rho_s, rho*v, E).
## Rows for momentum and energy are zero.
##
## Contract: jac is a caller-provided nvars*nvars scratch block which is
## zeroed and overwritten; callers targeting a larger matrix (e.g. the
## per-element blocks of a block-Jacobi preconditioner) must evaluate into
## a local block and add into the per-upt diagonal sub-block.  This is the
## Jacobian of the INSTANTANEOUS source (finite-rate, sub-steps = 0); it
## does not describe the sub-stepped source variants.  T is treated as
## exactly consistent with u; the fixed-count T_iter in the kernels makes
## this approximate, which is fine for preconditioning but means the
## result is not the exact Frechet derivative of the coded RHS.
<%pyfr:macro name='net_rate_jacobian' params='q, qh, jac'>
  for (int _i = 0; _i < ${nvars*nvars}; _i += 1)
    jac[_i] = 0.0;

  fpdtype_t rho = q[${rhoix}];
  fpdtype_t T = q[${Tix}];

  // Concentrations (log space)
  fpdtype_t log_cs[${ns}];
% for n in range(ns):
  log_cs[${n}] = log(fmax(${fpdtype_min}, rho*q[${n}]*${1.0/mcf[n].MW}));
% endfor

  // Gibbs free energy (non-dimensional G/RT) and its T derivative
  fpdtype_t gbs[${ns}];
  fpdtype_t dgbs[${ns}];
  fpdtype_t logT = log(T);
  fpdtype_t Tinv = 1.0/T;
  fpdtype_t log_prefRuT = ${math.log(101325.0/mcf.Ru)} + log(Tinv);
  fpdtype_t log_prefRuTinv = ${math.log(mcf.Ru/101325.0)} + logT;
% for n in range(ns):
  gbs[${n}] = ${mcf[n].gbs_expr('T', 'logT', 'Tinv')};
  dgbs[${n}] = -(${mcf[n].h_expr('T', False)})*Tinv*Tinv;
% endfor

  // Per-reaction d(rp)/dC accumulated into the species block;
  // sum(nu*d(rp)/dT) staged in the energy column
% for rtype in ('elementary', 'three-body', 'falloff-lindemann', 'falloff-troe'):
% for rxn in mcf.reactions_by_type(rtype):
  ${rxn.jac_block(nvars, vsmall=fpdtype_min)}
% endfor
% endfor

  // Assemble w.r.t. conserved variables via the temperature chain
  // dT/du = ((ke - e_j)/(rho cv), -v_d/(rho cv), 1/(rho cv))
  fpdtype_t _R = ${Rsum};
  fpdtype_t _rcv = 1.0/(rho*(qh[1] - _R));
  fpdtype_t _ke = 0.5*(${kesum});
  fpdtype_t _dTdr[${ns}];
% for n in range(ns):
  _dTdr[${n}] = (_ke - (qh[${4 + n}] - ${mcf.Ru/mcf[n].MW}*T))*_rcv;
% endfor
  fpdtype_t _MW[${ns}] = ${pyfr.carray(mws)};
  fpdtype_t _MWinv[${ns}] = ${pyfr.carray(mwinvs)};
  for (int _k = 0; _k < ${ns}; _k += 1)
  {
    fpdtype_t _dwT = jac[_k*${nvars} + ${Eix}];
    for (int _j = 0; _j < ${ns}; _j += 1)
      jac[_k*${nvars} + _j] = _MW[_k]*(jac[_k*${nvars} + _j]*_MWinv[_j]
                                       + _dwT*_dTdr[_j]);
% for i in range(ndims):
    jac[_k*${nvars} + ${vix + i}] = -_MW[_k]*_dwT*q[${vix + i}]*_rcv;
% endfor
    jac[_k*${nvars} + ${Eix}] = _MW[_k]*_dwT*_rcv;
  }
</%pyfr:macro>
