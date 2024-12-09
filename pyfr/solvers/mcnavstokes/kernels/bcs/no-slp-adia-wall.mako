<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor

% for i in range(ndims):
    ur[${i + vix}] = -ul[${i + vix}];
% endfor
    ur[${Eix}] = ul[${Eix}];

    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

% for n in range(ns):
    ur[${n}] = ul[${n}];
% endfor

% for i in range(ndims):
    ur[${i + vix}] = 0.0;
% endfor

    fpdtype_t rho = ${" + ".join([f"ul[{n}]" for n in range(ns)])};
    ur[${Eix}] = ul[${Eix}]
                     - (0.5/rho)*${pyfr.dot('ul[{i}]', i=(vix,vix + ndims))};

    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};

    // The LDG state is unique. We need to create an inconsistent state between
    // ur and qr. We clearly want to set the velocities to zero, and then compute
    // the primitive s based on just internal energy (no KE). But when we go to
    // compute T and species gradients => grad_ur, we use ul,grad_ul
    // in their full quantity, so when we go to viscous flux add, we need to use
    // the same values to get proper normal gradients on the wall. Further, on
    // a wall, we don't want to add the \tau*(ul-ur) term. Therefore, after
    // we compute the primitive s qr, we reset the conserved quantities, ur, to
    // be identical to ul.

% for i in range(ndims):
    ur[${i + vix}] = ul[${i + vix}];
% endfor
    ur[${Eix}] = ul[${Eix}];

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, ql, qhl, nl, grad_ul, grad_ur'>
    fpdtype_t rhoE = ul[${Eix}];
    fpdtype_t rho = ql[${rhoix}];
    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t E = rhoE*rcprho;

    // Copy all gradients to the right side, we will keep momentum, but we will
    // correct the energy and species terms such that the T and Y gradients
    // computed are orthogonal to the normal vector
% for i, j in pyfr.ndrange(ndims, nvars):
    grad_ur[${i}][${j}] = grad_ul[${i}][${j}];
% endfor

% if ndims == 2:

    fpdtype_t rho_x = ${" + ".join([f"grad_ul[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_ul[1][{n}]" for n in range(ns)])};

    // Velocity
    fpdtype_t u = ql[${vix + 0}];
    fpdtype_t v = ql[${vix + 1}];

    // Velocity derivatives (rho*d[u,v]/d[x,y])
    fpdtype_t u_x = grad_ul[0][${vix + 0}] - u*rho_x;
    fpdtype_t u_y = grad_ul[1][${vix + 0}] - u*rho_y;
    fpdtype_t v_x = grad_ul[0][${vix + 1}] - v*rho_x;
    fpdtype_t v_y = grad_ul[1][${vix + 1}] - v*rho_y;

    fpdtype_t rhoE_x = grad_ul[0][${Eix}];
    fpdtype_t rhoE_y = grad_ul[1][${Eix}];

    // Compute temperature derivatives (rho*cv*dT/d[x,y])
    fpdtype_t e_Y_Y_x;
    fpdtype_t e_Y_Y_y;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_y', 'ul', 'ql', 'qhl', 'grad_ul', 'rho_x', 'rho_y')};
    fpdtype_t T_x = rhoE_x - E*rho_x - u*u_x - v*v_x - rho*e_Y_Y_x;
    fpdtype_t T_y = rhoE_y - E*rho_y - u*u_y - v*v_y - rho*e_Y_Y_y;

    // Enforce no normal component of temperature gradient
    fpdtype_t Tdotn = T_x*nl[0] + T_y*nl[1];
    grad_ur[0][${Eix}] -= Tdotn*nl[0];
    grad_ur[1][${Eix}] -= Tdotn*nl[1];

    // Enforce zero normal species gradient in wall
    fpdtype_t Y_x, Y_y, Ydotn;
%   for n in range(ns):
    // Species derivative (rho*dY/d[x,y])
    Y_x = grad_ul[0][${n}] - ql[${n}]*rho_x;
    Y_y = grad_ul[1][${n}] - ql[${n}]*rho_y;
    Ydotn = Y_x*nl[0] + Y_y*nl[1];
    grad_ur[0][${n}] -= Ydotn*nl[0];
    grad_ur[1][${n}] -= Ydotn*nl[1];
%   endfor

% elif ndims == 3:
    fpdtype_t rho_x = ${" + ".join([f"grad_ul[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_ul[1][{n}]" for n in range(ns)])};
    fpdtype_t rho_z = ${" + ".join([f"grad_ul[2][{n}]" for n in range(ns)])};

    // Velocity
    fpdtype_t u = ql[${vix + 0}];
    fpdtype_t v = ql[${vix + 1}];
    fpdtype_t w = ql[${vix + 2}];

    // Velocity derivatives (rho*d[u,v,w]/d[x,y,z])
    fpdtype_t u_x = grad_ul[0][${vix + 0}] - u*rho_x;
    fpdtype_t u_y = grad_ul[1][${vix + 0}] - u*rho_y;
    fpdtype_t u_z = grad_ul[2][${vix + 0}] - u*rho_z;
    fpdtype_t v_x = grad_ul[0][${vix + 1}] - v*rho_x;
    fpdtype_t v_y = grad_ul[1][${vix + 1}] - v*rho_y;
    fpdtype_t v_z = grad_ul[2][${vix + 1}] - v*rho_z;
    fpdtype_t w_x = grad_ul[0][${vix + 2}] - w*rho_x;
    fpdtype_t w_y = grad_ul[1][${vix + 2}] - w*rho_y;
    fpdtype_t w_z = grad_ul[2][${vix + 2}] - w*rho_z;

    fpdtype_t rhoE_x = grad_ul[0][${Eix}];
    fpdtype_t rhoE_y = grad_ul[1][${Eix}];
    fpdtype_t rhoE_z = grad_ul[2][${Eix}];

    // Compute temperature derivatives (rho*cv*dT/d[x,y,z])
    fpdtype_t e_Y_Y_x;
    fpdtype_t e_Y_Y_y;
    fpdtype_t e_Y_Y_z;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_y', 'e_Y_Y_z', 'ul', 'ql', 'qhl', 'grad_ul', 'rho_x', 'rho_y', 'rho_z')};
    fpdtype_t T_x = rhoE_x - E*rho_x - u*u_x - v*v_x - w*w_x - rho*e_Y_Y_x;
    fpdtype_t T_y = rhoE_y - E*rho_y - u*u_y - v*v_y - w*w_y - rho*e_Y_Y_y;
    fpdtype_t T_z = rhoE_z - E*rho_z - u*u_z - v*v_z - w*w_z - rho*e_Y_Y_z;

    // Enforce no normal component of temperature gradient
    fpdtype_t Tdotn = T_x*nl[0] + T_y*nl[1] + T_z*nl[2];
    grad_ur[0][${Eix}] -= Tdotn*nl[0];
    grad_ur[1][${Eix}] -= Tdotn*nl[1];
    grad_ur[2][${Eix}] -= Tdotn*nl[2];

    // Enforce zero normal species gradient
    fpdtype_t Y_x, Y_y, Y_z, Ydotn;
%   for n in range(ns):
    // Species derivative (rho*dY/d[x,y,z])
    Y_x = grad_ul[0][${n}] - ql[${n}]*rho_x;
    Y_y = grad_ul[1][${n}] - ql[${n}]*rho_y;
    Y_z = grad_ul[2][${n}] - ql[${n}]*rho_z;
    Ydotn = Y_x*nl[0] + Y_y*nl[1] + Y_z*nl[2];
    grad_ur[0][${n}] -= Ydotn*nl[0];
    grad_ur[1][${n}] -= Ydotn*nl[1];
    grad_ur[2][${n}] -= Ydotn*nl[2];
%   endfor

% endif
</%pyfr:macro>
