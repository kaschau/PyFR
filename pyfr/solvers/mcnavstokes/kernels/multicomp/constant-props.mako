<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% ns, vix, Eix, rhoix, pix, Tix = pyfr.thermix(c['ns'], ndims) %>

<%pyfr:macro name='mixture_transport' params='u, q, qh, qt'>
  fpdtype_t rho = q[${rhoix}];

  // Mole fraction
  fpdtype_t X[${ns}];
  {
    fpdtype_t mass = 0.0;
% for n in range(ns):
    X[${n}] = q[${n}] * ${1.0 / c['MW'][n]};
    mass += X[${n}];
% endfor
    fpdtype_t invmass = 1.0 / mass;
% for n in range(ns):
    X[${n}] *= invmass;
% endfor
  }

  // Mixture viscosity
  fpdtype_t mu = 0.0;
% for n in range(ns):
    // ${c['names'][n]} viscosity
    {
      fpdtype_t phitemp = 0.0;
      fpdtype_t mu_n = ${c['mu0'][n]};
      % for n2 in range(ns):
        {
        fpdtype_t mu_n2 = ${max(c['mu0'][n2], fpdtype_eps)};
        fpdtype_t num = 1.0 + sqrt(mu_n / mu_n2 * sqrt(mu_n2 / mu_n));
        fpdtype_t phi = num*num * ${1.0/(math.sqrt(8.0) * math.sqrt(1.0 + c['MW'][n]/c['MW'][n2]))};
        phitemp += phi * X[${n2}];
        }
      % endfor
      mu += mu_n * X[${n}] / phitemp;
    }
% endfor

  qt[0] = mu;

  // Mixture thermal conductivity
  fpdtype_t kappa;
  {
    fpdtype_t sum1 = 0.0;
    fpdtype_t sum2 = 0.0;
    % for n in range(ns):
      // ${c['names'][n]} thermal conductivity
      sum1 += X[${n}] * ${c['kappa0'][n]};
      sum2 += X[${n}] * ${1.0 / (fpdtype_eps + c['kappa0'][n])};
    % endfor
    kappa = 0.5*(sum1 + 1.0 / sum2);
    qt[1] = kappa;
  }

  // Lewis number approximation
% for n in range(ns):
  qt[${2 + n}] = kappa/(rho * qh[1] * ${c['Le'][n]});
% endfor

</%pyfr:macro>