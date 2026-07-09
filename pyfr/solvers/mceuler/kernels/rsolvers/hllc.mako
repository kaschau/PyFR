<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
${pyfr.eos_check('hllc', fluid, 'mc-cpg', 'mc-tpg')}\
<%include file='pyfr.solvers.mceuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right primitive state
    ${fluid.decl('ul', 'rho, invrho, v, Y, p, a, gammamix', suffix='l')}
    ${fluid.decl('ur', 'rho, invrho, v, Y, p, a, gammamix', suffix='r')}

    // Compute the left and right fluxes
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t nf_fl, nf_fr, nf_fsl, nf_fsr;
    fpdtype_t va[${ndims}];
    fpdtype_t usl[${nvars}], usr[${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'pl', 'vl', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'pr', 'vr', 'fr')};

    fpdtype_t sqrtrl = sqrt(rhol);
    fpdtype_t sqrtrr = sqrt(rhor);

    // Get the normal left and right velocities
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // Compute the Roe-averaged velocity
    fpdtype_t nv = (sqrtrl*nvl + sqrtrr*nvr)/(sqrtrl + sqrtrr);

    // Compute the Roe-averaged enthalpy
    fpdtype_t H = (sqrtrl*(pr + ur[${nvars - 1}])
                 + sqrtrr*(pl + ul[${nvars - 1}]))
                / (sqrtrl*rhor + sqrtrr*rhol);

    fpdtype_t inv_rar = 1 / (sqrtrl + sqrtrr);
% for i in range(ndims):
    va[${i}] = (vl[${i}]*sqrtrl + vr[${i}]*sqrtrr)*inv_rar;
% endfor

    fpdtype_t qq = ${pyfr.dot('va[{i}]', i=ndims)};

    // Roe-averaged ratio of specific heats and speed of sound
    fpdtype_t gamma = (sqrtrl*gammamixl + sqrtrr*gammamixr)*inv_rar;
    fpdtype_t a = sqrt((gamma - 1)*(H - 0.5*qq));

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = min(nv - a, nvl - al);
    fpdtype_t sr = max(nv + a, nvr + ar);
    fpdtype_t sstar = (pr - pl + rhol*nvl*(sl - nvl)
                               - rhor*nvr*(sr - nvr)) /
                      (rhol*(sl - nvl) - rhor*(sr - nvr));

    // Star state common factors
    fpdtype_t ul_com = (sl - nvl) / (sl - sstar);
    fpdtype_t ur_com = (sr - nvr) / (sr - sstar);

    // Star state mass
    fpdtype_t rusl = ul_com*rhol;
    fpdtype_t rusr = ur_com*rhor;

    // Star state species
% for k in range(ns):
    usl[${k}] = rusl*Yl[${k}];
    usr[${k}] = rusr*Yr[${k}];
% endfor

    // Star state momentum
% for i in range(ndims):
    usl[${ns + i}] = rusl*(vl[${i}] + (sstar - nvl)*n[${i}]);
    usr[${ns + i}] = rusr*(vr[${i}] + (sstar - nvr)*n[${i}]);
%endfor

    // Star state energy
    usl[${nvars - 1}] = ul_com*(ul[${nvars - 1}] + (sstar - nvl)*
                                (rhol*sstar + pl/(sl - nvl)));
    usr[${nvars - 1}] = ur_com*(ur[${nvars - 1}] + (sstar - nvr)*
                                (rhor*sstar + pr/(sr - nvr)));

    // Output
% for i in range(nvars):
    nf_fl = ${' + '.join(f'n[{j}]*fl[{j}][{i}]' for j in range(ndims))};
    nf_fr = ${' + '.join(f'n[{j}]*fr[{j}][{i}]' for j in range(ndims))};
    nf_fsl = nf_fl + sl*(usl[${i}] - ul[${i}]);
    nf_fsr = nf_fr + sr*(usr[${i}] - ur[${i}]);
    nf[${i}] = (0 <= sl) ? nf_fl : (sl <= 0 && 0 <= sstar) ? nf_fsl :
               (sstar <= 0 && 0 <= sr) ? nf_fsr : nf_fr;
% endfor
</%pyfr:macro>
