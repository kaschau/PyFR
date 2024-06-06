<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.eos.cpg.p_from_rho_E'/>

<%pyfr:macro name='inviscid_flux' params='s, f, p, v'>
    fpdtype_t invrho = 1.0/s[0], rhoE = s[${ndims + 1}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Compute the pressure
    ${pyfr.expand('p_from_rho_E', 'p', 'rhoE', 'invrho', 'rhov')};

    // Density and energy fluxes
% for i in range(ndims):
    f[${i}][0] = rhov[${i}];
    f[${i}][${ndims + 1}] = (rhoE + p)*v[${i}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = rhov[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor

   // Species fluxes
% for i, k in pyfr.ndrange(ndims,ns-1):
    f[${i}][${k+ndims+2}] = v[${i}]*s[${k+ndims+2}];
% endfor

</%pyfr:macro>

<%pyfr:macro name='inviscid_flux_1d' params='s, f, p, v'>
    fpdtype_t invrho = 1.0/s[0], rhoE = s[${nvars - 1}];

    // Compute the velocities
% for i in range(ndims):
    v[${i}] = invrho*s[${i + 1}];
% endfor

    // Compute the pressure
    p = ${c['gamma'] - 1}*(rhoE - 0.5*invrho*${pyfr.dot('s[{i}]', i=(1, ndims + 1))});

    // Density and energy fluxes
    f[0] = s[1];
    f[${nvars - 1}] = (rhoE + p)*v[0];

    // Momentum fluxes
    f[1] = s[1]*v[0] + p;
% for j in range(1, ndims):
    f[${j + 1}] = s[1]*v[${j}];
% endfor

   // Species fluxes
% for k in range(ns-1):
    f[${k+3}] = s[${k+3}]*v[0];
% endfor

</%pyfr:macro>
