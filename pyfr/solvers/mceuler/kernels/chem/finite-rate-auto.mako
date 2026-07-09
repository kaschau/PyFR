<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.chem.net-rate-of-production'/>

<%pyfr:macro name='finite_rate_auto' params='t, u, ploc, src'>
    ${fluid.decl('u', 'rho, invrho, Y, T', suffix='ch')}

    fpdtype_t tmpSrc[${ns}];
% for n in range(ns):
    src[${n}] = 0.0;
% endfor

    fpdtype_t tRemaining = 1.0;
    for (int iter = 0; iter < ${max_subs} && tRemaining > ${fpdtype_eps};
         iter++)
    {
        ${pyfr.expand('net_rate_of_production', 'Ych', 'Tch', 'rhoch',
                      'tmpSrc')};

        // Adaptive sub-step size
        fpdtype_t tSubRatio = tRemaining;
% for n in range(ns):
% if participates[n]:
        {
            fpdtype_t s_ = copysign(1.0, tmpSrc[${n}]);
            fpdtype_t headroom_ = 0.5*((1.0 - s_)*Ych[${n}]
                                       + (1.0 + s_)*(1.0 - Ych[${n}]));
            fpdtype_t dtSubMax_ = rhoch*fmax(headroom_, ${fpdtype_eps})
                                / fmax(fabs(tmpSrc[${n}]), ${fpdtype_min});
            tSubRatio = fmin(tSubRatio, dtSubMax_*${1.0/dt});
        }
% endif
% endfor

        // Mixture cp at the current sub-step state
        fpdtype_t cp_ = 0.0;
% for n in range(ns):
        cp_ += (${fluid.species[n].cp_expr('Tch')})*Ych[${n}];
% endfor

        fpdtype_t tSub = tSubRatio*${dt};

        // Temperature and species sub-step
        fpdtype_t dTdt_ = 0.0;
% for n in range(ns):
% if participates[n]:
        {
            dTdt_ -= (${fluid.species[n].h_expr('Tch')})*tmpSrc[${n}];
            Ych[${n}] += tmpSrc[${n}]*invrhoch*tSub;
        }
% endif
        src[${n}] += tmpSrc[${n}]*tSubRatio;
% endfor

        dTdt_ /= cp_*rhoch;
        Tch += dTdt_*tSub;

        tRemaining -= tSubRatio;
    }

% for i in range(ndims):
    src[${ns + i}] = 0.0;
% endfor
    src[${nvars - 1}] = 0.0;
</%pyfr:macro>
