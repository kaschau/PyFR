<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='mixture_state' params='s, p, T, Y, Rmix, cpmix'>
    fpdtype_t rho = s[0];
    fpdtype_t invrho = 1.0/s[0];
    fpdtype_t rhoE = s[${ndims + 1}];

    Y[${ns-1}] = 1.0;
% for i in range(ns-1):
    Y[${i}] = s[${i+ndims+2}]*invrho;
    Y[${ns-1}] -= Y[${i}];
% endfor

    // Compute mixture properties
    Rmix = 0.0;
    cpmix = 0.0;
% for k in range(ns):
    Rmix += Y[${k}]*${props['MWinv'][k]};
    cpmix += Y[${k}]*${props['cp0'][k]};
% endfor
    Rmix *= ${props['Ru']};

    // Internal energy (per mass)
    fpdtype_t e = (rhoE - 0.5*invrho*${pyfr.dot('s[{i}]', i=(1,ndims+1))})*invrho;

    // Equilibrium p, T
    T = e / (cpmix - Rmix);
    p = rho * Rmix * T;

</%pyfr:macro>