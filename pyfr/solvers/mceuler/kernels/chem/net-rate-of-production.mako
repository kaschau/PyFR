<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%!
import math
%>

<%pyfr:macro name='net_rate_of_production' params='Y, T, rho, omega'>
% for n in range(ns):
    omega[${n}] = 0.0;
% endfor

    // Concentrations (log space)
    fpdtype_t log_cs[${ns}];
% for n in range(ns):
    log_cs[${n}] = log(fmax(${fpdtype_min},
                            rho*Y[${n}]*${1.0/fluid.species[n].MW}));
% endfor

    // Gibbs free energy (non-dimensional G/RT)
    fpdtype_t gbs[${ns}];
    fpdtype_t logT = log(T);
    fpdtype_t Tinv = 1.0/T;
    fpdtype_t log_prefRuT = ${math.log(101325.0/RU)} + log(Tinv);
    fpdtype_t log_prefRuTinv = ${math.log(RU/101325.0)} + logT;
% for n in range(ns):
    gbs[${n}] = ${fluid.species[n].gbs_expr('T', 'logT', 'Tinv')};
% endfor

    // Rate constants and rates of progress
% for rtype in ('elementary', 'three-body', 'falloff-lindemann', 'falloff-troe'):
% for rxn in fluid.reactions_by_type(rtype):
    ${rxn.rate_block()}
% endfor
% endfor

    // Convert to mass production rates
% for n in range(ns):
    omega[${n}] *= ${fluid.species[n].MW};
% endfor
</%pyfr:macro>
