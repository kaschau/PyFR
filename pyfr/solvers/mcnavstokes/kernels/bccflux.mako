<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<%include file='pyfr.solvers.mceuler.kernels.multicomp.${mcf.eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.${bctype}'/>

% if bccfluxstate:
<%include file='pyfr.solvers.mcnavstokes.kernels.bcs.${bccfluxstate}'/>
% endif

<% vix, Eix, rhoix, pix, Tix = mcf.mcix(ndims) %>\

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              artviscl='in view fpdtype_t'
              nl='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};


    // Compute left thermodynamic quantities
    fpdtype_t ql[${nvars + 2}];
    fpdtype_t qhl[${4 + mcf.ns}];
    ${pyfr.expand('stateFrom-cons', 'ul', 'ql', 'qhl')};

    ${pyfr.expand('bc_common_flux_state', 'ul', 'ql', 'qhl', 'gradul', 'artviscl', 'norm_nl', 'mag_nl')};
</%pyfr:kernel>
