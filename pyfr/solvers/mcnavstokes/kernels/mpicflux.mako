<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${eos}.e_Y_Y_x'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.multicomp.${trans}'/>
<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.mceuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.flux'/>

<% beta, tau = c['ldg-beta'], c['ldg-tau'] %>
<% ns = c['ns'] %>

<%pyfr:kernel name='mpicflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout mpi fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              gradur='in mpi fpdtype_t[${str(ndims)}][${str(nvars)}]'
              artviscl='in view fpdtype_t'
              artviscr='in mpi fpdtype_t'
              nl='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    // Compute left thermodynamic quantities
    fpdtype_t ql[${nvars + 1}];
    fpdtype_t qhl[${3 + ns}];
    ${pyfr.expand('stateFrom-cons', 'ul', 'ql', 'qhl')};

    // Compute right thermodynamic quantities
    fpdtype_t qr[${nvars + 1}];
    fpdtype_t qhr[${3 + ns}];
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}], fvcomm;
    ${pyfr.expand('rsolve', 'ul', 'ur', 'ql', 'qr', 'qhl', 'qhr', 'norm_nl', 'ficomm')};

% if beta != -0.5:
    fpdtype_t fvl[${ndims}][${nvars}] = {{0}};
    // Compute transport properties
    fpdtype_t qtl[${ns+2}];
    ${pyfr.expand('mixture_transport', 'ul', 'ql', 'qhl', 'qtl')};
    ${pyfr.expand('viscous_flux_add', 'ul', 'gradul', 'ql', 'qhl', 'qtl', 'fvl')};
    ${pyfr.expand('artificial_viscosity_add', 'gradul', 'fvl', 'artviscl')};
% endif

% if beta != 0.5:
    fpdtype_t fvr[${ndims}][${nvars}] = {{0}};
    // Compute transport properties
    fpdtype_t qtr[${ns+2}];
    ${pyfr.expand('mixture_transport', 'ur', 'qr', 'qhr', 'qtr')};
    ${pyfr.expand('viscous_flux_add', 'ur', 'gradur', 'qr', 'qhr', 'qtr', 'fvr')};
    ${pyfr.expand('artificial_viscosity_add', 'gradur', 'fvr', 'artviscr')};
% endif

% for i in range(nvars):
% if beta == -0.5:
    fvcomm = ${' + '.join(f'norm_nl[{j}]*fvr[{j}][{i}]' for j in range(ndims))};
% elif beta == 0.5:
    fvcomm = ${' + '.join(f'norm_nl[{j}]*fvl[{j}][{i}]' for j in range(ndims))};
% else:
    fvcomm = ${0.5 + beta}*(${' + '.join(f'norm_nl[{j}]*fvl[{j}][{i}]'
                                         for j in range(ndims))})
           + ${0.5 - beta}*(${' + '.join(f'norm_nl[{j}]*fvr[{j}][{i}]'
                                         for j in range(ndims))});
% endif
% if tau != 0.0:
    fvcomm += ${tau}*(ul[${i}] - ur[${i}]);
% endif

    ul[${i}] = mag_nl*(ficomm[${i}] + fvcomm);
% endfor
</%pyfr:kernel>
