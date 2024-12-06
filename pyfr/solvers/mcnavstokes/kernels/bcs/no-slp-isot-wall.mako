<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='bc_rsolve_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    // Set right primatives
    qr[${pix}] = ql[${pix}];
    // Set wall velocity
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + vix}] = -ql[${i + vix}] + 2*${c[v]};
% endfor
    // Set Temperature
    qr[${Tix}] = ${c['T']};

    // Set species
% for n in range(ns):
    qr[${n}] = ql[${n}];
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, ql, qhl, nl, ur, qr, qhr' externs='ploc, t'>

    // Set right primatives
    qr[${pix}] = ql[${pix}];
    // Set wall velocity
% for i, v in enumerate('uvw'[:ndims]):
    qr[${i + vix}] = -ql[${i + vix}] + 2*${c[v]};
% endfor
    // Set Temperature
    qr[${Tix}] = ${c['T']};

    // Set species
% for n in range(ns):
    qr[${n}] = ql[${n}];
% endfor

    ${pyfr.expand('stateFrom-prims', 'ur', 'qr', 'qhr')};

</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, ql, qhl, nl, grad_ul, grad_ur'>

    // Copy all gradients to the right side, we will keep momentum, but we will
    // species terms such that Y gradients
    // computed are orthogonal to the normal vector
% for i, j in pyfr.ndrange(ndims, nvars):
    grad_ur[${i}][${j}] = grad_ul[${i}][${j}];
% endfor

% if ndims == 2:
    fpdtype_t rho_x = ${" + ".join([f"grad_ul[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_ul[1][{n}]" for n in range(ns)])};

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

    // Enforce zero normal species gradient in wall
    fpdtype_t rho_x = ${" + ".join([f"grad_ul[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_ul[1][{n}]" for n in range(ns)])};
    fpdtype_t rho_z = ${" + ".join([f"grad_ul[2][{n}]" for n in range(ns)])};

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