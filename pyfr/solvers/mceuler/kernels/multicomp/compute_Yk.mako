<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='compute_Yk' params='Y, s, invrho'>
   Y[${ns-1}] = 1.0;
   fpdtype_t test_sum = 0.0;
% for i in range(ns-1):
   Y[${i}] = s[${i+ndims+2}]*invrho;
   Y[${ns-1}] -= Y[${i}];
   test_sum += Y[${i}];
% endfor
##TODO: Renormalize
</%pyfr:macro>