<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\
<% div = [1.0, 2.0, 3.0, 4.0] %>\

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    fpdtype_t rho = u[0];
    fpdtype_t T = q[${ndims + 1}];
    fpdtype_t lnT = log(T);
    e = 0.0;
    fpdtype_t Ymin = ${fpdtype_max};
    fpdtype_t Ysum = 0.0;
    // Compute mixture entropy
% for n in range(ns):
    // ${c['names'][n]} Entropy
    {
      fpdtype_t es;
      if (T < ${N7[n,0]})
      {
      <% m = 8 %>
          es = ${N7[n,m]*Ru/MW[n]}*lnT;
          es += T*(${'+ T*('.join(str(c) for c in N7[n,m+1:m+5]*Ru/MW[n]/div)+')'*3});
          es += ${N7[n,m+6]*Ru/MW[n]};
      }else
      {
      <% m = 1 %>
          es = ${N7[n,m]*Ru/MW[n]}*lnT;
          es += T*(${'+ T*('.join(str(c) for c in N7[n,m+1:m+5]*Ru/MW[n]/div)+')'*3});
          es += ${N7[n,m+6]*Ru/MW[n]};
      }
      e += es * q[${Yix+n}];
    }
%endfor
    e *= u[0];
    e = ((T > 0.0)) ? exp(e) : ${fpdtype_max};

</%pyfr:macro>