<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='mixture_transport' params='u, q, qh, qt'>
<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\
<% MW = c['MW'] %>\
<% muPoly = c['muPoly'] %>\
<% kappaPoly = c['kappaPoly'] %>\
## Dij is a symmetric matrix, so just store half and index appropriately
<% DijPoly = c['DijPoly'] %>\
<% Dijix = lambda n,n2 : int(ns * (ns-1) / 2 - (ns - n) * (ns - n - 1)/2 + n2) %>\
  fpdtype_t p = q[0];
  fpdtype_t T = q[${ndims+1}];

  fpdtype_t mu_sp[${ns}] = {0};
  fpdtype_t kappa_sp[${ns}] = {0};
  fpdtype_t invDij[${int((ns + 1)*ns/2)}] = {0};

  // Mole fraction
  fpdtype_t MWmix = 0.0;
  fpdtype_t X[${ns}];
  {
    fpdtype_t mass = 0.0;
% for n in range(ns):
    X[${n}] = q[${Yix+n}]*${1.0/MW[n]};
    mass += X[${n}];
% endfor
    fpdtype_t invmass = 1.0/mass;
% for n in range(ns):
    X[${n}] *= invmass;
    MWmix += X[${n}] * ${MW[n]};
% endfor
  }

  // Evaluate all property poly'l
  fpdtype_t logT = log(T);
  fpdtype_t sqrtT = sqrt(T);
  fpdtype_t sqrtsqrtT = sqrt(sqrtT);
  fpdtype_t T_3o2 = T*sqrtT;

% for n in range(ns):
// ${c['names'][n]} viscosity, thermal conductivity, diffusion coefficients
  mu_sp[${n}] = ${'+ logT*('.join(str(c) for c in muPoly[n,:])+')'*4};
  kappa_sp[${n}] = ${'+ logT*('.join(str(c) for c in kappaPoly[n,:])+')'*4};
  // Set to correct dimensions
  mu_sp[${n}] *= sqrtsqrtT;
  mu_sp[${n}] *= mu_sp[${n}];
  kappa_sp[${n}] *= sqrtT;
% for n2 in range(n, ns):
<% ix = Dijix(n,n2)%>\
    invDij[${ix}] = 1.0/((${'+ logT*('.join(str(c) for c in DijPoly[ix,:])+')'*4})*T_3o2);
% endfor
% endfor

  // Now every species' property is computed, generate mixture values
  // Mixture viscosity
  {
  fpdtype_t mu = 0.0;
% for n in range(ns):
    // ${c['names'][n]} viscosity
    {
    fpdtype_t phitemp = 0.0;
% for n2 in range(ns):
      {
      fpdtype_t num = 1.0 + sqrt(mu_sp[${n}] / mu_sp[${n2}] * ${math.sqrt(MW[n2] / MW[n])});
      fpdtype_t phi = num*num*${1.0/(math.sqrt(8.0) * math.sqrt(1.0 + MW[n]/MW[n2]))};
      phitemp += phi * X[${n2}];
      }
% endfor
    mu += mu_sp[${n}] * X[${n}] / phitemp;
    }
% endfor
    qt[0] = mu;
  }

  // Mixture thermal conductivity
  {
    fpdtype_t sum1 = 0.0;
    fpdtype_t sum2 = 0.0;
% for n in range(ns):
      // ${c['names'][n]} thermal conductivity
      sum1 += X[${n}] * kappa_sp[${n}];
      sum2 += X[${n}] / kappa_sp[${n}];
% endfor
    fpdtype_t kappa = 0.5*(sum1 + 1.0 / sum2);
    qt[1] = kappa;
  }

  // Mixture species diffusion coefficient
% for n in range(ns):
  {
    fpdtype_t sum1 = 0.0;
    fpdtype_t sum2 = 0.0;
% for n2 in range(ns):
% if n != n2:
##Symmetric
<% ix = Dijix(n,n2) if n2>=n else Dijix(n2,n)%>\
      sum1 += X[${n2}] * invDij[${ix}];
      sum2 += X[${n2}] * ${MW[n2]} * invDij[${ix}];
% endif
% endfor
    // account for pressure
    sum1 *= p;
    sum2 *= p * X[${n}] / (MWmix - ${MW[n]} * X[${n}] + 1e-12);
    // HACK
    fpdtype_t temp = sum1 + sum2;
    if (fabs(temp) > 1e-10){
      qt[${2+n}] = 1.0 / temp;
    }else{
      qt[${2+n}] = 0.0;
    }
  }
% endfor

#ifdef DEBUG
  printf("*********************************\n");
  printf("TRANSPORT PROPERTIES\n");
  printf("INPUT STATE\n");
  printf("trans&rho = %.14f\n", u[0]);
  printf("trans&p = %.14f\n", q[0]);
  printf("trans&T = %.14f\n", q[${ndims+1}]);
% for n in range(ns):
  printf("trans&Y_${c['names'][n]} = %.14f\n", q[${Yix+n}]);
% endfor

  printf("\nCOMPUTED PROPERTIES\n");
  printf("trans&MWmix = %.14f\n", MWmix);
  printf("trans&mu = %.14f\n", qt[0]);
  printf("trans&kappa = %.14f\n", qt[1]);
% for n in range(ns):
  printf("trans&D_${c['names'][n]} = %.14f\n", qt[${2+n}]);
% endfor
  printf("*********************************\n");

#endif
</%pyfr:macro>