<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>
<%include file='pyfr.solvers.mpeuler.kernels.bcs.${bctype}'/>

<% ns = c['ns'] %>

<%pyfr:kernel name='bccjump' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              jumpl='out view fpdtype_t[3]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});

    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('stateFrom-cons', 'ul', 'ql', 'qhl')};

    // Write out the jumps
    jumpl[0] = 0;
    jumpl[1] = fabs(ql[0]);
    jumpl[2] = mag_nl;
</%pyfr:kernel>