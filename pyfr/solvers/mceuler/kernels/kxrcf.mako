<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mceuler.kernels.multicomp.${eos}.stateFrom-cons'/>

<%pyfr:kernel name='kxrcf' ndim='1'
              jump='in fpdtype_t[${str(nfpts)}][5]'
              sensor='out fpdtype_t[2]'
              mass='in broadcast fpdtype_t[1][${str(nfpts)}]'>
    // Integrate jump around element
    sensor[0] = ${pyfr.dot('mass[0][{k}]', 'jump[{k}][0]', k=nfpts)};
    sensor[1] = ${pyfr.dot('mass[0][{k}]', 'jump[{k}][2]', k=nfpts)};
    fpdtype_t s_normp = jump[0][1];
    fpdtype_t s_normd = jump[0][3];
% for i in range(1, nfpts):
    s_normp = fmax(s_normp, jump[${i}][1]);
    s_normd = fmax(s_normd, jump[${i}][3]);
% endfor

    fpdtype_t sarea = ${pyfr.dot('mass[0][{k}]', 'jump[{k}][4]', k=nfpts)};
% if ndims == 2:
    fpdtype_t h = sarea*${1/math.pi};
% elif ndims == 3:
    fpdtype_t h = sqrt(sarea*${1/math.pi});
% endif

    sensor[0] = fabs(sensor[0])/(pow(h, ${order + 1.0})*sarea*s_normp);
    sensor[1] = fabs(sensor[1])/(pow(h, ${order + 1.0})*sarea*s_normd);
</%pyfr:kernel>