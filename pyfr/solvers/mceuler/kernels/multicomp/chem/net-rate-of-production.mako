<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%import math %>

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='net_rate_of_production' params='Y, T, rho, omega'>

  % for n in range(mcf.ns):
    omega[${n}] = 0.0;
  % endfor

  // Concentrations (log space)
  fpdtype_t log_cs[${mcf.ns}];
  % for n in range(mcf.ns):
    log_cs[${n}] = log(fmax(${fpdtype_min}, rho*Y[${n}]*${1.0/mcf[n].MW}));
  % endfor

  // Gibbs free energy (non-dimensional G/RT)
  fpdtype_t gbs[${mcf.ns}];
  fpdtype_t logT = log(T);
  fpdtype_t Tinv = 1.0/T;
  fpdtype_t log_prefRuT = ${math.log(101325.0/mcf.Ru)} + log(Tinv);
  fpdtype_t log_prefRuTinv = ${math.log(mcf.Ru/101325.0)} + logT;
  % for n in range(mcf.ns):
    gbs[${n}] = ${mcf[n].gbs_expr('T', 'logT', 'Tinv')};
  % endfor

  // Rate constants and rates of progress
  // Elementary reactions
  % for rxn in mcf.reactions_by_type('elementary'):
  ${rxn.rate_block()}
  % endfor

  // Three-body reactions
  % for rxn in mcf.reactions_by_type('three-body'):
  ${rxn.rate_block()}
  % endfor

  // Lindemann falloff reactions
  % for rxn in mcf.reactions_by_type('falloff-lindemann'):
  ${rxn.rate_block()}
  % endfor

  // Troe falloff reactions
  % for rxn in mcf.reactions_by_type('falloff-troe'):
  ${rxn.rate_block()}
  % endfor

  // Convert to mass production rates
  % for n in range(mcf.ns):
    omega[${n}] *= ${mcf[n].MW};
  % endfor
</%pyfr:macro>
