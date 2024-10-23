<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\

% if ndims == 2:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, q, qh, qt, fout'>
    fpdtype_t rho = q[${rhoix}];
    fpdtype_t rhou = uin[${vix}], rhov = uin[${vix + 1}];
    fpdtype_t rhoE = uin[${Eix}];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = q[${vix}], v = q[${vix + 1}];
    fpdtype_t mu = qt[0];
    fpdtype_t kappa = qt[1];

    fpdtype_t rho_x = ${" + ".join([f"grad_uin[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_uin[1][{n}]" for n in range(ns)])};

    // Velocity derivatives (d[u,v]/d[x,y])
    fpdtype_t u_x = rcprho*(grad_uin[0][${vix}] - u*rho_x);
    fpdtype_t u_y = rcprho*(grad_uin[1][${vix}] - u*rho_y);
    fpdtype_t v_x = rcprho*(grad_uin[0][${vix + 1}] - v*rho_x);
    fpdtype_t v_y = rcprho*(grad_uin[1][${vix + 1}] - v*rho_y);

    fpdtype_t rhoE_x = grad_uin[0][${Eix}];
    fpdtype_t rhoE_y = grad_uin[1][${Eix}];

    // Compute temperature derivatives (dT/d[x,y])
    fpdtype_t rcpcv = qh[0]/qh[1];
    fpdtype_t E = rhoE*rcprho;
    fpdtype_t e_Y_Y_x = 0.0;
    fpdtype_t e_Y_Y_y = 0.0;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_y', 'uin', 'q', 'qh', 'grad_uin')};
    fpdtype_t T_x = rcpcv*(rcprho*(rhoE_x - E*rho_x) - u*u_x - v*v_x - e_Y_Y_x);
    fpdtype_t T_y = rcpcv*(rcprho*(rhoE_y - E*rho_y) - u*u_y - v*v_y - e_Y_Y_y);

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*mu*(u_x - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_yy = -2*mu*(v_y - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_xy = -mu*(v_x + u_y);

    fout[0][${vix    }] += t_xx;  fout[1][${vix    }] += t_xy;
    fout[0][${vix + 1}] += t_xy;  fout[1][${vix + 1}] += t_yy;

    // Thermal diffusion
    fout[0][${Eix}] += u*t_xx + v*t_xy + -kappa*T_x;
    fout[1][${Eix}] += u*t_xy + v*t_yy + -kappa*T_y;

    // Species diffusion
    fpdtype_t Y_x, Y_y;
    fpdtype_t Jx, Jy;
%   for n in range(ns):
    // Species derivative (rho*dY/d[x,y])
      Y_x = grad_uin[0][${n}] - q[${n}]*rho_x;
      Y_y = grad_uin[1][${n}] - q[${n}]*rho_y;

      // Species mass diffusion
      Jx = qt[${2 + n}]*Y_x;
      fout[0][${n}] -= Jx;
      Jy = qt[${2 + n}]*Y_y;
      fout[1][${n}] -= Jy;

      // Species thermal diffusion
      fout[0][${Eix}] -= qh[${4 + n}] * Jx;
      fout[1][${Eix}] -= qh[${4 + n}] * Jy;
%   endfor
</%pyfr:macro>

% elif ndims == 3:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, q, qh, qt, fout'>
    fpdtype_t rho = q[${rhoix}];
    fpdtype_t rhou = uin[${vix}], rhov = uin[${vix + 1}], rhow = uin[${vix + 2}];
    fpdtype_t rhoE = uin[${Eix}];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = q[${vix}], v = q[${vix + 1}], w = q[${vix + 2}];
    fpdtype_t mu = qt[0];
    fpdtype_t kappa = qt[1];

    fpdtype_t rho_x = ${" + ".join([f"grad_uin[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"grad_uin[1][{n}]" for n in range(ns)])};
    fpdtype_t rho_z = ${" + ".join([f"grad_uin[2][{n}]" for n in range(ns)])};

    // Velocity derivatives (grad[u,v,w])
    fpdtype_t u_x = rcprho*(grad_uin[0][${vix    }] - u*rho_x);
    fpdtype_t u_y = rcprho*(grad_uin[1][${vix    }] - u*rho_y);
    fpdtype_t u_z = rcprho*(grad_uin[2][${vix    }] - u*rho_z);
    fpdtype_t v_x = rcprho*(grad_uin[0][${vix + 1}] - v*rho_x);
    fpdtype_t v_y = rcprho*(grad_uin[1][${vix + 1}] - v*rho_y);
    fpdtype_t v_z = rcprho*(grad_uin[2][${vix + 1}] - v*rho_z);
    fpdtype_t w_x = rcprho*(grad_uin[0][${vix + 2}] - w*rho_x);
    fpdtype_t w_y = rcprho*(grad_uin[1][${vix + 2}] - w*rho_y);
    fpdtype_t w_z = rcprho*(grad_uin[2][${vix + 2}] - w*rho_z);

    fpdtype_t rhoE_x = grad_uin[0][${Eix}];
    fpdtype_t rhoE_y = grad_uin[1][${Eix}];
    fpdtype_t rhoE_z = grad_uin[2][${Eix}];

    // Compute temperature derivatives (dT/d[x,y,z])
    fpdtype_t rcpcv = qh[0]/qh[1];
    fpdtype_t E = rhoE*rcprho;
    fpdtype_t e_Y_Y_x = 0.0;
    fpdtype_t e_Y_Y_y = 0.0;
    fpdtype_t e_Y_Y_z = 0.0;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_y', 'e_Y_Y_z', 'uin', 'q', 'qh', 'grad_uin')};
    fpdtype_t T_x = rcpcv*(rcprho*(rhoE_x - E*rho_x) - u*u_x - v*v_x - w*w_x - e_Y_Y_x);
    fpdtype_t T_y = rcpcv*(rcprho*(rhoE_y - E*rho_y) - u*u_y - v*v_y - w*w_y - e_Y_Y_y);
    fpdtype_t T_z = rcpcv*(rcprho*(rhoE_z - E*rho_z) - u*u_z - v*v_z - w*w_z - e_Y_Y_z);

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*mu*(u_x - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_yy = -2*mu*(v_y - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_zz = -2*mu*(w_z - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_xy = -mu*(v_x + u_y);
    fpdtype_t t_xz = -mu*(u_z + w_x);
    fpdtype_t t_yz = -mu*(w_y + v_z);

    fout[0][${vix + 0}] += t_xx;     fout[1][${vix + 0}] += t_xy;     fout[2][${vix + 0}] += t_xz;
    fout[0][${vix + 1}] += t_xy;     fout[1][${vix + 1}] += t_yy;     fout[2][${vix + 1}] += t_yz;
    fout[0][${vix + 2}] += t_xz;     fout[1][${vix + 2}] += t_yz;     fout[2][${vix + 2}] += t_zz;

    fout[0][${Eix}] += u*t_xx + v*t_xy + w*t_xz + -kappa*T_x;
    fout[1][${Eix}] += u*t_xy + v*t_yy + w*t_yz + -kappa*T_y;
    fout[2][${Eix}] += u*t_xz + v*t_yz + w*t_zz + -kappa*T_z;

    // Species diffusion
    fpdtype_t Y_x, Y_y, Y_z;
    fpdtype_t Jx, Jy, Jz;
%   for n in range(ns):
      // Species derivative (rho*dY/d[x,y,z])
      Y_x = grad_uin[0][${n}] - q[${n}]*rho_x;
      Y_y = grad_uin[1][${n}] - q[${n}]*rho_y;
      Y_z = grad_uin[2][${n}] - q[${n}]*rho_z;

      // Species mass diffusion
      Jx = qt[${2 + n}]*Y_x;
      fout[0][${n}] -= Jx;
      Jy = qt[${2 + n}]*Y_y;
      fout[1][${n}] -= Jy;
      Jz = qt[${2 + n}]*Y_z;
      fout[2][${n}] -= Jz;

      // Species thermal diffusion
      fout[0][${Eix}] -= qh[${4 + n}] * Jx;
      fout[1][${Eix}] -= qh[${4 + n}] * Jy;
      fout[2][${Eix}] -= qh[${4 + n}] * Jz;
%   endfor
</%pyfr:macro>
% endif
