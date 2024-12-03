<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

% if ndims == 2:
<%pyfr:macro name='e_Y_Y_x' params='e_Y_Y_x, e_Y_Y_y, u, q, qh, gradu, rho_x, rho_y'>
    fpdtype_t invrho = 1.0/q[${rhoix}];
    fpdtype_t T = q[${Tix}];

    e_Y_Y_x = 0.0;
    e_Y_Y_y = 0.0;

% for n in range(ns):
    {
      fpdtype_t Y_x =  invrho*(gradu[0][${n}] - q[${n}]*rho_x);
      fpdtype_t Y_y =  invrho*(gradu[1][${n}] - q[${n}]*rho_y);
      fpdtype_t e_Y = T*(${c['cp0'][n] - c['Ru']/c['MW'][n]});
      e_Y_Y_x += e_Y * Y_x;
      e_Y_Y_y += e_Y * Y_y;
    }
% endfor

</%pyfr:macro>
% elif ndims == 3:
<%pyfr:macro name='e_Y_Y_x' params='e_Y_Y_x, e_Y_Y_y, e_Y_Y_z, u, q, qh, gradu, rho_x, rho_y, rho_z'>
    fpdtype_t invrho = 1.0/q[${rhoix}];
    fpdtype_t T = q[${Tix}];

    e_Y_Y_x = 0.0;
    e_Y_Y_y = 0.0;
    e_Y_Y_z = 0.0;
% for n in range(ns):
    {
      fpdtype_t Y_x =  invrho*(gradu[0][${n}] - q[${n}]*rho_x);
      fpdtype_t Y_y =  invrho*(gradu[1][${n}] - q[${n}]*rho_y);
      fpdtype_t Y_z =  invrho*(gradu[2][${n}] - q[${n}]*rho_z);
      fpdtype_t e_Y = T*(${c['cp0'][n] - c['Ru']/c['MW'][n]});
      e_Y_Y_x += e_Y * Y_x;
      e_Y_Y_y += e_Y * Y_y;
      e_Y_Y_z += e_Y * Y_z;
    }
% endfor
</%pyfr:macro>
% endif