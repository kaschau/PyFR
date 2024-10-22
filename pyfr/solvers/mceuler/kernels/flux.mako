<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='inviscid_flux' params='u, f, q'>
    fpdtype_t rho = q[${rhoix}];
    fpdtype_t p = q[${pix}];
    fpdtype_t invrho = 1.0/rho;
    fpdtype_t rhoE = u[${Eix}];

   // Species fluxes
% for i, n in pyfr.ndrange(ndims,ns):
   f[${i}][${n}] = q[${vix + i}]*u[${n}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${vix + j}] = u[${vix + i}]*q[${vix + j}]${' + p' if i == j else ''};
% endfor

    // Energy flux
% for i in range(ndims):
    f[${i}][${Eix}] = (rhoE + p)*q[${vix + i}];
% endfor

</%pyfr:macro>

<%pyfr:macro name='inviscid_flux_1d' params='u, f, q'>
    fpdtype_t invrho = 1.0/q[${rhoix}];
    fpdtype_t rhoE = u[${Eix}];
    fpdtype_t p = q[${pix}];

    // Compute the velocities
    fpdtype_t v[${ndims}];
% for i in range(ndims):
    v[${i}] = q[${vix + i}];
% endfor

    // Species fluxes
% for n in range(ns):
    f[${n}] = u[${n}]*v[0];
% endfor

    // Momentum fluxes
    f[${vix}] = u[${vix}]*v[0] + q[${pix}];
% for j in range(1, ndims):
    f[${j + vix}] = u[${vix}]*v[${j}];
% endfor

    // Energy fluxes
    f[${Eix}] = (rhoE + q[${pix}])*v[0];

</%pyfr:macro>
