<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>

<% Yix = ndims + 2 %>\
<% ns = c['ns'] %>\

<%pyfr:macro name='rsolve_1d' params='ul, ur, ql, qr, qhl, qhr, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr, fsl, fsr;
    fpdtype_t usl[${nvars}], usr[${nvars}];

    ${pyfr.expand('inviscid_flux_1d', 'ul', 'fl', 'ql')};
    ${pyfr.expand('inviscid_flux_1d', 'ur', 'fr', 'qr')};

    % for i in range(ndims):
      vl[${i}] = ql[${i+1}];
      vr[${i}] = qr[${i+1}];
    % endfor
    pl = ql[0];
    pr = qr[0];

    fpdtype_t sqrtrl = sqrt(ul[0]);
    fpdtype_t sqrtrr = sqrt(ur[0]);

    // Compute the Roe-averaged enthalpy
    fpdtype_t H = (sqrtrl*(pr + ur[${ndims + 1}])
                 + sqrtrr*(pl + ul[${ndims + 1}]))
                / (sqrtrl*ur[0] + sqrtrr*ul[0]);

    // Roe average sound speed
    fpdtype_t u = (sqrtrl*vl[0] + sqrtrr*vr[0]) /
                  (sqrtrl + sqrtrr);

    fpdtype_t gamma = (sqrtrl*qhl[0] + sqrtrr*qhr[0]) /
                      (sqrtrl + sqrtrr);
    fpdtype_t a = sqrt((gamma - 1)*(H - 0.5*u*u));

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = u - a;
    fpdtype_t sr = u + a;
    fpdtype_t sstar = (pr - pl + ul[0]*vl[0]*(sl - vl[0])
                               - ur[0]*vr[0]*(sr - vr[0])) /
                      (ul[0]*(sl - vl[0]) - ur[0]*(sr - vr[0]));

    // Star state common factors
    fpdtype_t ul_com = (sl - vl[0]) / (sl - sstar);
    fpdtype_t ur_com = (sr - vr[0]) / (sr - sstar);

    // Star state mass
    usl[0] = ul_com*ul[0];
    usr[0] = ur_com*ur[0];

    // Star state momenetum
    usl[1] = ul_com*ul[0]*sstar;
    usr[1] = ur_com*ur[0]*sstar;
% for i in range(2, ndims + 1):
    usl[${i}] = ul_com*ul[${i}];
    usr[${i}] = ur_com*ur[${i}];
% endfor

    // Star state energy
    usl[${ndims + 1}] = ul_com*(ul[${ndims + 1}] + (sstar - vl[0])*
                                (ul[0]*sstar + pl/(sl - vl[0])));
    usr[${ndims + 1}] = ur_com*(ur[${ndims + 1}] + (sstar - vr[0])*
                                (ur[0]*sstar + pr/(sr - vr[0])));

    // Star state species
% for n in range(ns - 1):
    usl[${Yix + n}] = usl[0]*ql[${Yix + n}];
    usr[${Yix + n}] = usr[0]*qr[${Yix + n}];
% endfor

    // Output
% for i in range(nvars):
    fsl = fl[${i}] + sl*(usl[${i}] - ul[${i}]);
    fsr = fr[${i}] + sr*(usr[${i}] - ur[${i}]);
    nf[${i}] = (0 <= sl) ? fl[${i}] : (sl <= 0 && 0 <= sstar) ? fsl :
               (sstar <= 0 && 0 <= sr) ? fsr : fr[${i}];
% endfor
</%pyfr:macro>

<%include file='pyfr.solvers.mceuler.kernels.rsolvers.rsolve1d'/>
