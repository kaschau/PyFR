<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${eos}.e_Y_Y_x'/>

% if ndims == 2:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, q, qh, qt, fout'>
    fpdtype_t rho = uin[0], rhou = uin[1], rhov = uin[2], rhoE = uin[3];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov;
    fpdtype_t mu = qt[0];
    fpdtype_t kappa = qt[1];

    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];

    // Velocity derivatives (d[u,v]/d[x,y])
    fpdtype_t u_x = rcprho*(grad_uin[0][1] - u*rho_x);
    fpdtype_t u_y = rcprho*(grad_uin[1][1] - u*rho_y);
    fpdtype_t v_x = rcprho*(grad_uin[0][2] - v*rho_x);
    fpdtype_t v_y = rcprho*(grad_uin[1][2] - v*rho_y);

    fpdtype_t rhoE_x = grad_uin[0][3];
    fpdtype_t rhoE_y = grad_uin[1][3];

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

    fout[0][1] += t_xx;     fout[1][1] += t_xy;
    fout[0][2] += t_xy;     fout[1][2] += t_yy;

    // Thermal diffusion
    fout[0][3] += u*t_xx + v*t_xy + -kappa*T_x;
    fout[1][3] += u*t_xy + v*t_yy + -kappa*T_y;

    // Species diffusion
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>
    fpdtype_t Y_x, Y_y;
    ## TODO diffusion correction term
%   for n in range(ns-1):
    // Species derivative (rho*dY/d[x,y])
    Y_x =  grad_uin[0][${Yix+n}] - q[${Yix+n}]*rho_x;
    Y_y =  grad_uin[1][${Yix+n}] - q[${Yix+n}]*rho_y;

    fout[0][${Yix+n}] += -qt[${2+n}] * Y_x;
    fout[1][${Yix+n}] += -qt[${2+n}] * Y_y;

    ## TODO species thermal diffusion
%   endfor

</%pyfr:macro>
% elif ndims == 3:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
    fpdtype_t rho  = uin[0];
    fpdtype_t rhou = uin[1], rhov = uin[2], rhow = uin[3];
    fpdtype_t rhoE    = uin[4];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov, w = rcprho*rhow;

    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];
    fpdtype_t rho_z = grad_uin[2][0];

    // Velocity derivatives (grad[u,v,w])
    fpdtype_t u_x = rcprho*(grad_uin[0][1] - u*rho_x);
    fpdtype_t u_y = rcprho*(grad_uin[1][1] - u*rho_y);
    fpdtype_t u_z = rcprho*(grad_uin[2][1] - u*rho_z);
    fpdtype_t v_x = rcprho*(grad_uin[0][2] - v*rho_x);
    fpdtype_t v_y = rcprho*(grad_uin[1][2] - v*rho_y);
    fpdtype_t v_z = rcprho*(grad_uin[2][2] - v*rho_z);
    fpdtype_t w_x = rcprho*(grad_uin[0][3] - w*rho_x);
    fpdtype_t w_y = rcprho*(grad_uin[1][3] - w*rho_y);
    fpdtype_t w_z = rcprho*(grad_uin[2][3] - w*rho_z);

    fpdtype_t rhoE_x = grad_uin[0][4];
    fpdtype_t rhoE_y = grad_uin[1][4];
    fpdtype_t rhoE_z = grad_uin[2][4];

    // Compute temperature derivatives (dT/d[x,y,z])
    fpdtype_t rcpcv = qh[0]/qh[1];
    fpdtype_t E = rhoE*rcprho;
    fpdtype_t e_Y_Y_x = 0.0;
    fpdtype_t e_Y_Y_y = 0.0;
    fpdtype_t e_Y_Y_z = 0.0;
    ${pyfr.expand('e_Y_Y_x', 'e_Y_Y_x', 'e_Y_Y_z', 'uin', 'q', 'qh', 'grad_uin')};
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

    fout[0][1] += t_xx;     fout[1][1] += t_xy;     fout[2][1] += t_xz;
    fout[0][2] += t_xy;     fout[1][2] += t_yy;     fout[2][2] += t_yz;
    fout[0][3] += t_xz;     fout[1][3] += t_yz;     fout[2][3] += t_zz;

    fout[0][4] += u*t_xx + v*t_xy + w*t_xz + -kappa*T_x;
    fout[1][4] += u*t_xy + v*t_yy + w*t_yz + -kappa*T_y;
    fout[2][4] += u*t_xz + v*t_yz + w*t_zz + -kappa*T_z;

    // Species diffusion
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>
    fpdtype_t Y_x, Y_y;
    ## TODO diffusion correction term
%   for n in range(ns-1):
    // Species derivative (rho*dY/d[x,y])
    Y_x =  grad_uin[0][${Yix+n}] - q[${Yix+n}]*rho_x;
    Y_y =  grad_uin[1][${Yix+n}] - q[${Yix+n}]*rho_y;
    Y_z =  grad_uin[2][${Yix+n}] - q[${Yix+n}]*rho_z;

    fout[0][${Yix+n}] += -qt[${2+n}] * Y_x;
    fout[1][${Yix+n}] += -qt[${2+n}] * Y_y;
    fout[2][${Yix+n}] += -qt[${2+n}] * Y_z;

    ## TODO species thermal diffusion
%   endfor
</%pyfr:macro>
% endif