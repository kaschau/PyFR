<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.chem.net-rate-of-production'/>

<% tSub = dt / float(sub_steps) %>

<%pyfr:macro name='finite_rate_substep' params='t, u, ploc, src'>
    ${fluid.decl('u', 'rho, invrho, Y, T', suffix='ch')}

    fpdtype_t tmpSrc[${ns}];
% for n in range(ns):
    src[${n}] = 0.0;
% endfor

    for (int nSub = 0; nSub < ${sub_steps}; nSub++)
    {
        ${pyfr.expand('net_rate_of_production', 'Ych', 'Tch', 'rhoch',
                      'tmpSrc')};

        // Mixture cp at the current sub-step state
        fpdtype_t cp_ = 0.0;
% for n in range(ns):
        cp_ += (${fluid.species[n].cp_expr('Tch')})*Ych[${n}];
% endfor

        // Temperature sub-step
        fpdtype_t dTdt_ = 0.0;
% for n in range(ns):
% if participates[n]:
        dTdt_ -= (${fluid.species[n].h_expr('Tch')})*tmpSrc[${n}];
% endif
        src[${n}] += tmpSrc[${n}]*${tSub/dt};
% endfor
        dTdt_ /= cp_*rhoch;
        Tch += dTdt_*${tSub};

        // Species sub-step (clamp + normalise)
        fpdtype_t Ysum_ = 0.0;
% for n in range(ns):
% if participates[n]:
        Ych[${n}] = fmin(1.0, fmax(0.0, Ych[${n}]
                                        + tmpSrc[${n}]*invrhoch*${tSub}));
% endif
        Ysum_ += Ych[${n}];
% endfor
        fpdtype_t Ysuminv_ = 1.0/Ysum_;
% for n in range(ns):
        Ych[${n}] *= Ysuminv_;
% endfor
    }

% for i in range(ndims):
    src[${ns + i}] = 0.0;
% endfor
    src[${nvars - 1}] = 0.0;
</%pyfr:macro>
