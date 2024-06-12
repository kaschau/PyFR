<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.baseadvec.kernels.smats'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.mixture_state'/>

<% smats = 'smats_l' if 'linear' in ktype else 'smats' %>

<%pyfr:kernel name='tflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              f='out fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              verts='in broadcast-col fpdtype_t[${str(nverts)}][${str(ndims)}]'
              upts='in broadcast-row fpdtype_t[${str(ndims)}]'>
% if 'linear' in ktype:
    // Compute the S matrices
    fpdtype_t ${smats}[${ndims}][${ndims}], djac;
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', smats, 'djac')};
% endif

    // Compute thermodynamic properties
    fpdtype_t q[${nvars+1}];
    fpdtype_t qh[${3}];
    ${pyfr.expand('mixture_state', 'u', 'q', 'qh')};

    // Compute the flux
    fpdtype_t ftemp[${ndims}][${nvars}];
    fpdtype_t v[${ndims}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp', 'q')};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join(f'{smats}[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
