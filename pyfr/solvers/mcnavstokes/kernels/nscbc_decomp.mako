<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% invsq2 = 2**-0.5 %>
<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='WU_dot_div-cartesian' params='div, Phi, q, qh'>
  fpdtype_t rho = q[${rhoix}];
  fpdtype_t invrho = 1.0/rho;
  fpdtype_t c = qh[2];
  fpdtype_t invc = 1.0/c;
  fpdtype_t invcsq = invc*invc;

  fpdtype_t nx = norm_nl[0];
  fpdtype_t ny = norm_nl[1];
  fpdtype_t v[${ndims}] = {${','.join([f'q[{vix+i}]' for i in range(ndims)])}};

  fpdtype_t gmo = qh[0] - 1.0;
  fpdtype_t k = 0.5*(${pyfr.dot('v[{i}]', i=ndims)});

  fpdtype_t divr = ${'+'.join(f'div[{n}]' for n in range(mcf.ns))};
  fpdtype_t invRmix = 1.0/(qh[1] - qh[1]/qh[0]);
  fpdtype_t sumdhdivYi = 0.0;
  % for n in range(mcf.ns):
  {
    fpdtype_t dhi = qh[${4+n}] - invRmix*qh[1]*q[${Tix}]*${mcf.Ru/mcf[n].MW};
    sumdhdivYi += div[${n}]*dhi;
  }
  % endfor

  % for n in range(mcf.ns):
    Phi[${ndims+n}] = invrho*(div[${n}] - divr*q[${n}]);
  % endfor

% if ndims == 2:

  Phi[0] = divr - gmo*invcsq*(divr*k - div[${vix}]*v[0] - div[${vix+1}]*v[1] + div[${Eix}] - sumdhdivYi);
  Phi[1] = invrho*(divr*(ny*v[0] - nx*v[1]) - div[${vix}]*ny + div[${vix+1}]*nx);
  Phi[${nvars-1}] = ${invsq2}*invc*invrho*(gmo*(div[${Eix}] - sumdhdivYi + divr*k - (div[${vix}]*v[0] + div[${vix+1}]*v[1])) + c*(div[${vix}]*nx + div[${vix+1}]*ny - divr*(nx*v[0] + ny*v[1])));
  Phi[${nvars  }] = ${invsq2}*invc*invrho*(gmo*(div[${Eix}] - sumdhdivYi + divr*k - (div[${vix}]*v[0] + div[${vix+1}]*v[1])) - c*(div[${vix}]*nx + div[${vix+1}]*ny - divr*(nx*v[0] + ny*v[1])));

% elif ndims == 3:

  fpdtype_t nz = norm_nl[2];

  fpdtype_t temp = -divr*k + div[${vix}]*v[0] + div[${vix+1}]*v[1] + div[${vix+2}]*v[2] - div[${Eix}] + sumdhdivYi;
  Phi[0] = gmo*nx*invcsq*(temp) + invrho*(-div[${vix+2}]*ny + div[${vix+1}]*nz + divr*(nx*rho - nz*v[1] + ny*v[2]));
  Phi[1] = gmo*ny*invcsq*(temp) + invrho*( div[${vix+2}]*nx - div[${vix  }]*nz + divr*(ny*rho + nz*v[0] - nx*v[2]));
  Phi[2] = gmo*nz*invcsq*(temp) + invrho*(-div[${vix+1}]*nx + div[${vix  }]*ny + divr*(nz*rho - ny*v[0] + nx*v[1]));
  Phi[${nvars-1}] = invc*invrho*${invsq2}*(gmo*(div[${Eix}] + divr*k - sumdhdivYi - (div[${vix}]*v[0] + div[${vix+1}]*v[1] + div[${vix+2}]*v[2])) + c*(div[${vix}]*nx + div[${vix+1}]*ny + div[${vix+2}]*nz - divr*(nx*v[0] + ny*v[1] + nz*v[2])));
  Phi[${nvars  }] = invc*invrho*${invsq2}*(gmo*(div[${Eix}] + divr*k - sumdhdivYi - (div[${vix}]*v[0] + div[${vix+1}]*v[1] + div[${vix+2}]*v[2])) - c*(div[${vix}]*nx + div[${vix+1}]*ny + div[${vix+2}]*nz - divr*(nx*v[0] + ny*v[1] + nz*v[2])));

% endif
</%pyfr:macro>


<%pyfr:macro name='WUinv_dot_Phi-cartesian' params='Phi, div, q, qh'>
  fpdtype_t rho = q[${rhoix}];
  fpdtype_t invrho = 1.0/rho;
  fpdtype_t c = qh[2];
  fpdtype_t invc = 1.0/c;

  fpdtype_t nx = norm_nl[0];
  fpdtype_t ny = norm_nl[1];

  fpdtype_t gmo = qh[0] - 1.0;

  fpdtype_t v[${ndims}] = {${','.join([f'q[{vix+i}]' for i in range(ndims)])}};
  fpdtype_t k = 0.5*(${pyfr.dot('v[{i}]', i=ndims)});

  fpdtype_t term1 = 0.0;
  fpdtype_t invRmix = 1.0/(qh[1] - qh[1]/qh[0]);
  % for n in range(mcf.ns):
  {
    fpdtype_t dh = qh[${4+n}] - invRmix*qh[1]*q[${Tix}]*${mcf.Ru/mcf[n].MW};
    term1 += Phi[${ndims+n}]*dh;
  }
  % endfor
  fpdtype_t term2 = k + qh[3]/rho - qh[1]/qh[0]*q[${Tix}];

  % if ndims == 2:

  fpdtype_t temp = Phi[0] + invc*${invsq2}*rho*(Phi[${nvars-1}] + Phi[${nvars}]);
  % for n in range(mcf.ns):
    div[${n}] = Phi[${ndims+n}]*rho + temp*q[${n}];
  % endfor

  div[${vix+0}] = -Phi[1]*ny*rho + Phi[0]*v[0] + ${invsq2}*rho*invc*(Phi[${nvars-1}]*(c*nx + v[0]) - Phi[${nvars}]*(c*nx - v[0]));
  div[${vix+1}] =  Phi[1]*nx*rho + Phi[0]*v[1] + ${invsq2}*rho*invc*(Phi[${nvars-1}]*(c*ny + v[1]) - Phi[${nvars}]*(c*ny - v[1]));
  div[${Eix}] = rho*term1 + Phi[0]*term2
               - Phi[1]*rho*(ny*v[0] - nx*v[1]) +
               ${invsq2}*rho*(
                              (Phi[${nvars-1}] + Phi[${nvars}])*(c/gmo + term2*invc) +
                              (Phi[${nvars-1}] - Phi[${nvars}])*(nx*v[0] + ny*v[1])
                             );

  % elif ndims == 3:
  fpdtype_t nz = norm_nl[2];

  fpdtype_t temp = Phi[0]*nx + Phi[1]*ny + Phi[2]*nz + ${invsq2}*rho*invc*(Phi[${nvars-1}] + Phi[${nvars}]);
  % for n in range(mcf.ns):
    div[${n}] = Phi[${ndims+n}]*rho + temp*q[${n}];
  % endfor

  div[${vix+0}] = rho*(Phi[2]*ny - Phi[1]*nz) + Phi[0]*nx*v[0] + Phi[1]*ny*v[0] + Phi[2]*nz*v[0] + ${invsq2}*invc*rho*(Phi[${nvars-1}]*(c*nx+v[0]) + Phi[${nvars}]*(-c*nx + v[0]));
  div[${vix+1}] = rho*(Phi[0]*nz - Phi[2]*nx) + Phi[0]*nx*v[1] + Phi[1]*ny*v[1] + Phi[2]*nz*v[1] + ${invsq2}*invc*rho*(Phi[${nvars-1}]*(c*ny+v[1]) + Phi[${nvars}]*(-c*ny + v[1]));
  div[${vix+2}] = rho*(Phi[1]*nx - Phi[0]*ny) + Phi[0]*nx*v[2] + Phi[1]*ny*v[2] + Phi[2]*nz*v[2] + ${invsq2}*invc*rho*(Phi[${nvars-1}]*(c*nz+v[2]) + Phi[${nvars}]*(-c*nz + v[2]));
  div[${Eix  }] = Phi[0]*(term2*nx + nz*rho*v[1] - ny*rho*v[2]) +
                  Phi[1]*(term2*ny - nz*rho*v[0] + nx*rho*v[2]) +
                  Phi[2]*(term2*nz + ny*rho*v[0] - nx*rho*v[1]) +
                rho*invc*${invsq2}/gmo*(Phi[${nvars-1}]*(c*c + gmo*term2 + c*gmo*(nx*v[0] + ny*v[1] + nz*v[2]))
                                      + Phi[${nvars  }]*(c*c + gmo*term2 - c*gmo*(nx*v[0] + ny*v[1] + nz*v[2])));

  % endif
</%pyfr:macro>


<%pyfr:macro name='WU_dot_div-normal' params='div, Phi, q, qh'>
  fpdtype_t rho = q[${rhoix}];
  fpdtype_t invrho = 1.0/rho;
  fpdtype_t c = qh[2];
  fpdtype_t invc = 1.0/c;
  fpdtype_t invcsq = invc*invc;
  fpdtype_t gmo = qh[0] - 1.0;
  fpdtype_t v[${ndims}] = {${','.join([f'q[{vix+i}]' for i in range(ndims)])}};
  fpdtype_t k = 0.5*(${pyfr.dot('v[{i}]', i=ndims)});

  fpdtype_t un = ${pyfr.dot('norm_nl[{i}]','v[{i}]', i=ndims)};
  fpdtype_t ut1 = ${pyfr.dot('t1[{i}]','v[{i}]', i=ndims)};
  % if ndims == 3:
  fpdtype_t ut2 = ${pyfr.dot('t2[{i}]','v[{i}]', i=ndims)};
  % endif

  fpdtype_t div_v = ${'+'.join([f'div[{vix+i}]*v[{i}]' for i in range(ndims)])};
  fpdtype_t div_norm = ${'+'.join([f'div[{vix+i}]*norm_nl[{i}]' for i in range(ndims)])};
  fpdtype_t div_t1 = ${'+'.join([f'div[{vix+i}]*t1[{i}]' for i in range(ndims)])};
  % if ndims == 3:
  fpdtype_t div_t2 = ${'+'.join([f'div[{vix+i}]*t2[{i}]' for i in range(ndims)])};
  % endif

  fpdtype_t invRmix = 1.0/(qh[1] - qh[1]/qh[0]);
  fpdtype_t divr = ${'+'.join(f'div[{n}]' for n in range(mcf.ns))};
  fpdtype_t sumdhdivYi = 0.0;
  % for n in range(mcf.ns):
  {
    fpdtype_t dhi = qh[${4+n}] - invRmix*qh[1]*q[${Tix}]*${mcf.Ru/mcf[n].MW};
    sumdhdivYi += div[${n}]*dhi;
  }
  % endfor

  Phi[0] = divr - gmo*invcsq*(divr*k - div_v + div[${Eix}] - sumdhdivYi);
  % for i in range(ndims-1):
    Phi[${i + 1}] = invrho*(div_t${i+1} - divr*ut${i+1});
  % endfor
  % for n in range(mcf.ns):
  Phi[${ndims+n}] = invrho*(div[${n}] - divr*q[${n}]);
  % endfor
  Phi[${nvars-1}] = ${invsq2}*invc*invrho*(gmo*(divr*k - div_v + div[${Eix}] - sumdhdivYi) + c*(div_norm - divr*un));
  Phi[${nvars  }] = ${invsq2}*invc*invrho*(gmo*(divr*k - div_v + div[${Eix}] - sumdhdivYi) - c*(div_norm - divr*un));
</%pyfr:macro>


<%pyfr:macro name='WUinv_dot_Phi-normal' params='Phi, div, q, qh'>
  fpdtype_t rho = q[${rhoix}];
  fpdtype_t invrho = 1.0/rho;
  fpdtype_t c = qh[2];
  fpdtype_t invc = 1.0/c;
  fpdtype_t gmo = qh[0] - 1.0;
  fpdtype_t v[${ndims}] = {${','.join([f'q[{vix+i}]' for i in range(ndims)])}};
  fpdtype_t k = 0.5*(${pyfr.dot('v[{i}]', i=ndims)});

  fpdtype_t un = ${pyfr.dot('norm_nl[{i}]','v[{i}]', i=ndims)};
  fpdtype_t ut1 = ${pyfr.dot('t1[{i}]','v[{i}]', i=ndims)};
  % if ndims == 3:
  fpdtype_t ut2 = ${pyfr.dot('t2[{i}]','v[{i}]', i=ndims)};
  % endif

  fpdtype_t term1 = 0.0;
  fpdtype_t invRmix = 1.0/(qh[1] - qh[1]/qh[0]);
  % for n in range(mcf.ns):
  {
    fpdtype_t dh = qh[${4+n}] - invRmix*qh[1]*q[${Tix}]*${mcf.Ru/mcf[n].MW};
    term1 += Phi[${ndims+n}]*dh;
  }
  % endfor
  fpdtype_t term2 = k + qh[3]/rho - qh[1]/qh[0]*q[${Tix}];

  fpdtype_t temp = Phi[0] + invc*${invsq2}*rho*(Phi[${nvars-1}] + Phi[${nvars}]);
  % for n in range(mcf.ns):
    div[${n}] = Phi[${ndims+n}]*rho + temp*q[${n}];
  % endfor
  % for i in range(ndims):
    div[${vix+i}] = Phi[0]*v[${i}] + rho*(${'+'.join([f'Phi[{j+1}]*t{j+1}[{i}]' for j in range(ndims-1)])}) + ${invsq2}*rho*invc*(Phi[${nvars-1}]*(c*norm_nl[${i}] + v[${i}]) - Phi[${nvars}]*(c*norm_nl[${i}] - v[${i}]));
  % endfor

  div[${Eix}] = term2*Phi[0] +
               rho*(${'+'.join([f'Phi[{j+1}]*ut{j+1}' for j in range(ndims-1)])}) +
               rho*term1 +
               ${invsq2}*rho*(
                              (Phi[${nvars-1}] + Phi[${nvars}])*(c/gmo + term2*invc) +
                              (Phi[${nvars-1}] - Phi[${nvars}])*un
                             );
</%pyfr:macro>


<%pyfr:macro name='WUinv_dot_Phi-noop' params='Phi, div, q, qh'>
% for var in range(nvars):
div[${var}] = Phi[${var}];
% endfor
</%pyfr:macro>
<%pyfr:macro name='WU_dot_div-noop' params='div, Phi, q, qh'>
% for var in range(nvars):
Phi[${var}] = div[${var}];
% endfor
</%pyfr:macro>
