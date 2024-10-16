<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>

% if ndims == 2:
<%pyfr:macro name='e_Y_Y_x' params='e_Y_Y_x, e_Y_Y_y, u, q, qh, gradu'>
    fpdtype_t invrho = 1.0/u[0];
    fpdtype_t T = q[${ndims+1}];
    fpdtype_t rho_x = gradu[0][0];
    fpdtype_t rho_y = gradu[1][0];

    fpdtype_t Yns_x = 0.0;
    fpdtype_t Yns_y = 0.0;
% for n in range(ns-1):
    {
      fpdtype_t Y_x =  invrho*(gradu[0][${Yix+n}] - q[${Yix+n}]*rho_x);
      fpdtype_t Y_y =  invrho*(gradu[1][${Yix+n}] - q[${Yix+n}]*rho_y);
      fpdtype_t e_Y = T*(${c['cp0'][n] - c['Ru']/c['MW'][n]});
      e_Y_Y_x += e_Y * Y_x;
      e_Y_Y_y += e_Y * Y_y;
      Yns_x -= Y_x;
      Yns_y -= Y_y;
    }
% endfor
    // Compute Yns component
    {
      fpdtype_t e_Y = T*(${c['cp0'][ns-1] - c['Ru']/c['MW'][ns-1]});
      e_Y_Y_x += e_Y * Yns_x;
      e_Y_Y_y += e_Y * Yns_y;
    }

</%pyfr:macro>
% elif ndims == 3:
<%pyfr:macro name='e_Y_Y_x' params='e_Y_Y_x, e_Y_Y_y, e_Y_Y_z, u, q, qh, gradu'>
    fpdtype_t invrho = 1.0/u[0];
    fpdtype_t T = q[${ndims+1}];
    fpdtype_t rho_x = gradu[0][0];
    fpdtype_t rho_y = gradu[1][0];
    fpdtype_t rho_z = gradu[2][0];

    fpdtype_t Yns_x = 0.0;
    fpdtype_t Yns_y = 0.0;
    fpdtype_t Yns_z = 0.0;
% for n in range(ns-1):
    {
      fpdtype_t Y_x =  invrho*(gradu[0][${Yix+n}] - q[${Yix+n}]*rho_x);
      fpdtype_t Y_y =  invrho*(gradu[1][${Yix+n}] - q[${Yix+n}]*rho_y);
      fpdtype_t Y_z =  invrho*(gradu[2][${Yix+n}] - q[${Yix+n}]*rho_z);
      fpdtype_t e_Y = T*(${c['cp0'][n] - c['Ru']/c['MW'][n]});
      e_Y_Y_x += e_Y * Y_x;
      e_Y_Y_y += e_Y * Y_y;
      e_Y_Y_z += e_Y * Y_z;
      Yns_x -= Y_x;
      Yns_y -= Y_y;
      Yns_z -= Y_z;
    }
% endfor
    // Compute Yns component
    {
      fpdtype_t e_Y = T*(${c['cp0'][ns-1] - c['Ru']/c['MW'][ns-1]});
      e_Y_Y_x += e_Y * Yns_x;
      e_Y_Y_y += e_Y * Yns_y;
      e_Y_Y_z += e_Y * Yns_z;
    }
</%pyfr:macro>
% endif