<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Isotropic linear blend toward the cell mean:  u(alpha) = u + alpha*(uavg - u).
//
// alpha = 0 is the identity; alpha = 1 collapses every upt to the cell mean.
// Cell mean is preserved by construction (no mass-defect correction needed).
//
// Two variants: _density only touches u[*][0]; _full touches all variables.

<%pyfr:macro name='iso_lin_blend_density' params='u, uavg, alpha'>
    % for uidx in range(nupts):
    u[${uidx}][0] += alpha*(uavg[0] - u[${uidx}][0]);
    % endfor
</%pyfr:macro>

<%pyfr:macro name='iso_lin_blend_full' params='u, uavg, alpha'>
    % for uidx, vidx in pyfr.ndrange(nupts, nvars):
    u[${uidx}][${vidx}] += alpha*(uavg[${vidx}] - u[${uidx}][${vidx}]);
    % endfor
</%pyfr:macro>
