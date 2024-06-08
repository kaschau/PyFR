<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='p_from_rho_E' params='p, E, invrho, rhov'>
    // Internal energy (per mass)
    fpdtype_t e = (E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)})*invrho;

    // Compute mixture properties
    fpdtype_t Rmix = 0.0;
    fpdtype_t cp = 0.0;
% for k in range(ns):
    Rmix += Y[${k}]*${props['MWinv'][k]};
    cp += Y[${k}]*${props['cp0'][k]};
% endfor
    Rmix *= ${props['Ru']};

    p = s[0] * Rmix * e / (cp - Rmix);

</%pyfr:macro>
