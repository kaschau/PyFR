<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<% ns = c['ns'] %>

<%pyfr:kernel name='mpicjump' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              ur='in mpi fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              jumpl='out view fpdtype_t[5]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});

    fpdtype_t ql[${nvars+1}];
    fpdtype_t qhl[${3+ns}];
    ${pyfr.expand('stateFrom-cons', 'ul', 'ql', 'qhl')};

    fpdtype_t qr[${nvars+1}];
    fpdtype_t qhr[${3+ns}];
    ${pyfr.expand('stateFrom-cons', 'ur', 'qr', 'qhr')};

    // Write out the jumps
    jumpl[0] = mag_nl*(ql[0] - qr[0]);
    jumpl[1] = fabs(ql[0]);
    jumpl[2] = mag_nl*(ul[0] - ur[0]);
    jumpl[3] = fabs(ul[0]);
    jumpl[4] = mag_nl;
</%pyfr:kernel>