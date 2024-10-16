<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.baseadvec.kernels.smats'/>
<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${eos}.e_Y_Y_x'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${trans}'/>
<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.baseadvecdiff.kernels.transform_grad'/>
<%include file='pyfr.solvers.mceuler.kernels.flux'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.flux'/>

<% gradu = 'gradu' if 'fused' in ktype else 'f' %>
<% smats = 'smats_l' if 'linear' in ktype else 'smats' %>
<% rcpdjac = 'rcpdjac_l' if 'linear' in ktype else 'rcpdjac' %>

<%pyfr:kernel name='tflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              artvisc='in broadcast-col fpdtype_t'
              f='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'
              gradu='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t'
              verts='in broadcast-col fpdtype_t[${str(nverts)}][${str(ndims)}]'
              upts='in broadcast-row fpdtype_t[${str(ndims)}]'>
% if 'linear' in ktype:
    // Compute the S matrices
    fpdtype_t ${smats}[${ndims}][${ndims}], djac;
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', smats, 'djac')};
    fpdtype_t ${rcpdjac} = 1 / djac;
% endif
<% ns = c['ns'] %>

% if 'fused' in ktype:
    // Transform the corrected gradient
    ${pyfr.expand('transform_grad', gradu, smats, rcpdjac)};
% endif

    // Compute the flux (F = Fi + Fv)

    // Compute thermodynamic properties
    fpdtype_t q[${nvars + 1}];
    fpdtype_t qh[${4 + ns}];
    ${pyfr.expand('stateFrom-cons', 'u', 'q', 'qh')};

    fpdtype_t ftemp[${ndims}][${nvars}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp', 'q')};

    // Compute transport properties
    fpdtype_t qt[${2 + ns}];
    ${pyfr.expand('mixture_transport', 'u', 'q', 'qh', 'qt')};
    ${pyfr.expand('viscous_flux_add', 'u', gradu, 'q', 'qh', 'qt', 'ftemp')};
    ${pyfr.expand('artificial_viscosity_add', gradu, 'ftemp', 'artvisc')};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join(f'{smats}[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
