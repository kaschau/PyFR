<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>\

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% div = [1.0, 2.0, 3.0, 4.0, 5.0] %>\


% if ndims == 2:
<%pyfr:macro name='e_Y_Y_x' params='e_Y_Y_x, e_Y_Y_y, u, q, qh, gradu'>
    fpdtype_t invrho = 1.0/q[${rhoix}];
    fpdtype_t T = q[${Tix}];
    fpdtype_t rho_x = ${" + ".join([f"gradu[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"gradu[1][{n}]" for n in range(ns)])};

% for n in range(ns):
    {
      fpdtype_t Y_x =  invrho*(gradu[0][${n}] - q[${n}]*rho_x);
      fpdtype_t Y_y =  invrho*(gradu[1][${n}] - q[${n}]*rho_y);
      fpdtype_t hs = qh[${4 + n}];
      fpdtype_t e_Y = hs - T*(${c['Ru']/c['MW'][n]});
      e_Y_Y_x += e_Y * Y_x;
      e_Y_Y_y += e_Y * Y_y;
    }
% endfor

</%pyfr:macro>
% elif ndims == 3:
<%pyfr:macro name='e_Y_Y_x' params='e_Y_Y_x, e_Y_Y_y, e_Y_Y_z, u, q, qh, gradu'>
    fpdtype_t invrho = 1.0/q[${rhoix}];
    fpdtype_t T = q[${Tix}];
    fpdtype_t rho_x = ${" + ".join([f"gradu[0][{n}]" for n in range(ns)])};
    fpdtype_t rho_y = ${" + ".join([f"gradu[1][{n}]" for n in range(ns)])};
    fpdtype_t rho_z = ${" + ".join([f"gradu[2][{n}]" for n in range(ns)])};

% for n in range(ns):
    {
      fpdtype_t Y_x =  invrho*(gradu[0][${n}] - q[${n}]*rho_x);
      fpdtype_t Y_y =  invrho*(gradu[1][${n}] - q[${n}]*rho_y);
      fpdtype_t Y_z =  invrho*(gradu[2][${n}] - q[${n}]*rho_z);
      fpdtype_t hs = qh[${4 + n}];
      fpdtype_t e_Y = hs - T*(${c['Ru']/c['MW'][n]});
      e_Y_Y_x += e_Y * Y_x;
      e_Y_Y_y += e_Y * Y_y;
      e_Y_Y_z += e_Y * Y_z;
    }
% endfor

</%pyfr:macro>
% endif