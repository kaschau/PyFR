<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='mixture_transport' params='u, q, qh, qt'>
<% Yix = ndims + 2 %>
<% ns = c['ns'] %>
  fpdtype_t p = q[0];
  fpdtype_t T = q[${ndims+1}];

  fpdtype_t mu_sp[${ns}] = {0};
  fpdtype_t kappa_sp[${ns}] = {0};
  fpdtype_t Dij[${ns}][${ns}] = {0};

  // Mole fraction
  fpdtype_t MWmix = 0.0;
  fpdtype_t X[${ns}];
  {
    fpdtype_t mass = 0.0;
% for n in range(ns):
    X[${n}] = q[${Yix+n}]*${1.0/c['MW'][n]};
    mass += X[${n}];
% endfor
    fpdtype_t invmass = 1.0/mass;
% for n in range(ns):
    X[${n}] *= invmass;
    MWmix += X[${n}] * ${c['MW'][n]};
% endfor
  }

  // Evaluate all property poly'l
<% deg = 4 %>
  fpdtype_t logT = log(T);
  fpdtype_t sqrtT = exp(0.5 * logT);
  fpdtype_t T_3o2 = pow(T,1.5);
  fpdtype_t logT_n[${deg+1}];
  logT_n[0] = 1.0;
% for i in range(1,deg+1):
   logT_n[${i}] = logT * logT_n[${i-1}];
% endfor

<% muPoly = c['muPoly'] %>
<% kappaPoly = c['kappaPoly'] %>
<% DijPoly = c['DijPoly'] %>
% for n in range(ns):
// ${c['names'][n]} visc, kappa, Dij
% for i in range(deg+1):
  mu_sp[${n}] += ${muPoly[n, i]} * logT_n[${i}];
  kappa_sp[${n}] += ${kappaPoly[n, i]} * logT_n[${i}];
% for n2 in range(n, ns):
<% indx = int(ns * (ns-1) / 2 - (ns - n) * (ns - n - 1)/2 + n2) %>
  Dij[${n}][${n2}] += ${DijPoly[indx, i]} * logT_n[${i}];
% endfor
% endfor
  // Set to correct dimensions
  mu_sp[${n}] *= sqrtT;
  kappa_sp[${n}] *= sqrtT;
% for n2 in range(n, ns):
  Dij[${n}][${n2}] *= T_3o2;
% if n != n2:
  Dij[${n2}][${n}] = Dij[${n}][${n2}];
% endif

% endfor
% endfor

  // Now every species' property is computed, generate mixture values

  // Mixture viscosity
  {
  fpdtype_t mu = 0.0;
  fpdtype_t phi;
  fpdtype_t phitemp;
% for n in range(ns):
    // ${c['names'][n]} viscosity
    phitemp = 0.0;
% for n2 in range(ns):
      phi = pow((1.0 + sqrt(mu_sp[${n}] / mu_sp[${n2}] * sqrt(mu_sp[${n2}] / mu_sp[${n}]))), 2.0)
                    / (sqrt(8.0) * sqrt(1.0 + ${c['MW'][n]}/${c['MW'][n2]}));
      phitemp += phi * X[${n2}];
% endfor
    mu += mu_sp[${n}] * X[${n}] / phitemp;
% endfor
    qt[0] = mu;
  }

  // Mixture thermal conductivity
  fpdtype_t kappa;
  {
    fpdtype_t sum1 = 0.0;
    fpdtype_t sum2 = 0.0;
% for n in range(ns):
      // ${c['names'][n]} thermal conductivity
      sum1 += X[${n}] * kappa_sp[${n}];
      sum2 += X[${n}] / kappa_sp[${n}];
% endfor
    kappa = 0.5*(sum1 + 1.0 / sum2);
    qt[1] = kappa;
  }

  // mixture species diffusion coefficient
  {
    fpdtype_t sum1;
    fpdtype_t sum2;
% for n in range(ns):
    sum1 = 0.0;
    sum2 = 0.0;
% for n2 in range(ns):
% if n != n2:
      sum1 += X[${n2}] / Dij[${n}][${n2}];
      sum2 += X[${n2}] * ${c['MW'][n2]} / Dij[${n}][${n2}];
% endif
% endfor
  // Account for pressure
  sum1 *= p;
  sum2 *= p * X[${n}] / (MWmix - ${c['MW'][n]} * X[${n}]);
  qt[${2+n}] = 1.0 / (sum1 + sum2);
% endfor
  }
</%pyfr:macro>