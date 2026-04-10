<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>

<%pyfr:macro name='mixture_transport' params='u, q, qh, qt'>
  fpdtype_t rho = q[${rhoix}];

  // Mole fraction
  fpdtype_t X[${mcf.ns}];
  {
    fpdtype_t mass = 0.0;
% for n in range(mcf.ns):
    X[${n}] = q[${n}] * ${1.0 / mcf[n].MW};
    mass += X[${n}];
% endfor
    fpdtype_t invmass = 1.0 / mass;
% for n in range(mcf.ns):
    X[${n}] *= invmass;
% endfor
  }

  // Mixture viscosity
  fpdtype_t mu = 0.0;
% for n in range(mcf.ns):
    {
      fpdtype_t phitemp = 0.0;
      fpdtype_t mu_n = ${mcf[n].mu0};
      % for n2 in range(mcf.ns):
        {
        fpdtype_t mu_n2 = ${max(mcf[n2].mu0, fpdtype_eps)};
        fpdtype_t num = 1.0 + sqrt(mu_n / mu_n2 * sqrt(mu_n2 / mu_n));
        fpdtype_t phi = num*num * ${1.0/(math.sqrt(8.0) * math.sqrt(1.0 + mcf[n].MW/mcf[n2].MW))};
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
    % for n in range(mcf.ns):
      sum1 += X[${n}] * ${mcf[n].kappa0};
      sum2 += X[${n}] * ${1.0 / (fpdtype_eps + mcf[n].kappa0)};
    % endfor
    kappa = 0.5*(sum1 + 1.0 / sum2);
    qt[1] = kappa;
  }

  // Lewis number approximation
% for n in range(mcf.ns):
  qt[${2 + n}] = kappa/(rho * qh[1] * ${mcf[n].Le});
% endfor

</%pyfr:macro>
