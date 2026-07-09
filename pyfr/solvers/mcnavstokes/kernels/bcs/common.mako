<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_common_grad_zero' params='ub, nl, grad_ul, grad_ur'>
% for i, j in pyfr.ndrange(ndims, nvars):
    grad_ur[${i}][${j}] = 0;
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_grad_copy' params='ub, nl, grad_ul, grad_ur'>
% for i, j in pyfr.ndrange(ndims, nvars):
    grad_ur[${i}][${j}] = grad_ul[${i}][${j}];
% endfor
</%pyfr:macro>
