<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<% N7 = c['NASA7'] %>\
<% Ru = c['Ru'] %>\
<% MW = c['MW'] %>\
<% div = [1.0, 2.0, 3.0, 4.0] %>\

<%pyfr:macro name='compute_entropy' params='u, q, e'>

    fpdtype_t rho = q[${rhoix}];
    fpdtype_t T = q[${Tix}];
    fpdtype_t lnT = log(T);
    e = 0.0;
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
      e += es * q[${n}];
    }

    // Return the specific thermodynamic entropy (mass basis)
    e = (T > 0) ? e : ${fpdtype_max};
% endfor
</%pyfr:macro>