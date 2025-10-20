<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.backends.${_backend.name}.ikp'/>

<% se0 = math.log10(c['s0']/order**4) %>

## IKP version: kernel developer writes per-element logic
## Backend automatically transforms to block-level with cache blocking
<%pyfr:ikpkernel name='shocksensor_autoblock' ndim='1'
                 u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
                 artvisc='out fpdtype_t'>

    // Declare local arrays (per-element view)
    fpdtype_t ucol[${nupts}];
    int temp = 1;

    // Extract column (rho should be encoded as 1.0 + x + 10*y + 100*z)
    ${pyfr.expand('loadv', 'ucol', 'u', indices=[(i, svar) for i in range(nupts)])}

    fpdtype_t modes[${nupts}];
    ${pyfr.ikpexpand('gemv', 'modes', 'ucol', A=invvdm)}

    // Compute sensor
    fpdtype_t totEn = 0.0, pnEn = 1e-15;

% for i, bmode in enumerate(ind_modes):
    {
        fpdtype_t tmp = modes[${i}];
        totEn += tmp*tmp;
%   if bmode:
        pnEn += tmp*tmp;
%   endif
    }
% endfor

    fpdtype_t se  = ${1/math.log(10)}*log(pnEn/totEn);
    temp += 1;

    // Compute cell-wise artificial viscosity
    fpdtype_t mu = (se < ${se0 - c['kappa']})
                 ? 0.0
                 : ${0.5*c['max-artvisc']}*(1.0 + sin(${0.5*math.pi/c['kappa']}*(se - ${se0})));
    mu = (se < ${se0 + c['kappa']}) ? mu : ${c['max-artvisc']};


    artvisc = mu;
    printf("artvisc %e\n", artvisc);
    printf("temp %d\n", temp);
</%pyfr:ikpkernel>
