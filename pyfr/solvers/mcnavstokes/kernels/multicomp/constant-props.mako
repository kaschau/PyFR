<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='mixture_transport' params='u, q, qh, qt'>
  fpdtype_t rho = u[0];

  // Mole fraction
<% Yix = ndims + 2 %>
  fpdtype_t X[${ns}];
  {
    fpdtype_t mass = 0.0;
% for n in range(ns):
    X[${n}] = q[${Yix+n}]*${1.0/props['MW'][n]};
    mass += X[${n}];
% endfor
% for n in range(ns):
    X[${n}] /= mass;
% endfor
  }

  // Mixture viscosity
  {
  fpdtype_t mu = 0.0;
% for n in range(ns):
    fpdtype_t phitemp = 0.0;
    fpdtype_t mu_n = ${props['mu0'][n]};
    fpdtype_t MW_n = ${props['MW'][n]};
% for n2 in range(ns):
      fpdtype_t mu_n2 = ${props['mu0'][n2]};
      fpdtype_t MW_n2 = ${props['MW'][n2]};
      fpdtype_t phi = pow((1.0 + sqrt(mu_n / mu_n2 * sqrt(mu_n2 / mu_n))), 2.0)
                     / (sqrt(8.0) * sqrt(1.0 + MW_n/MW_n2));
      phitemp += phi * X[${n2}];
% endfor
    mu += mu_n * X[${n}] / phitemp;
% endfor
    qt[0] = mu;
  }

  // Mixture viscosity
  fpdtype_t kappa = 0.0;
  {
    fpdtype_t sum1 = 0.0;
    fpdtype_t sum2 = 0.0;
% for n in range(ns):
      sum1 += X[${ns}] * ${props['kappa0'][n]};
      sum2 += X[${ns}] / ${props['kappa0'][n]};
% endfor
    kappa = 0.5*(sum1 + 1.0/ sum2);
    qt[2] = kappa;
  }

  // Unity Lewis approximation
% for n in range(ns):
  qt[${2+n}] = kappa/(rho*qh[1]*${props['Le'][n]});
% endfor

</%pyfr:macro>