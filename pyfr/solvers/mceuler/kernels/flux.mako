<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% ns = c['ns'] %>
<% Yix = ndims + 2 %>

<%pyfr:macro name='inviscid_flux' params='u, f, q'>
    fpdtype_t rho = u[0];
    fpdtype_t p = q[0];
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${ndims + 1}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}], v[${ndims}];
% for i in range(ndims):
    rhov[${i}] = u[${i + 1}];
    v[${i}] = q[${i+1}];
% endfor

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
% for i, n in pyfr.ndrange(ndims,ns-1):
   f[${i}][${Yix+n}] = v[${i}]*u[${Yix+n}];
% endfor

</%pyfr:macro>

<%pyfr:macro name='inviscid_flux_1d' params='u, f, q'>
    fpdtype_t rho = u[0];
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${ndims + 1}];

    // Compute the velocities
    fpdtype_t v[${ndims}];
% for i in range(ndims):
    v[${i}] = q[${i + 1}];
% endfor

    // Density and energy fluxes
    f[0] = u[1];
    f[${ndims + 1}] = (rhoE + q[0])*v[0];

    // Momentum fluxes
    f[1] = u[1]*v[0] + q[0];
% for j in range(1, ndims):
    f[${j + 1}] = u[1]*v[${j}];
% endfor

    // Species fluxes
% for k in range(ns-1):
    f[${Yix+k}] = u[${Yix+k}]*v[0];
% endfor

</%pyfr:macro>
