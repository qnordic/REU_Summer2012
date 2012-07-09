function ImportGenRTR
   GenRTRBase=pwd;
   fprintf('Adding GenRTR paths from %s...\n',GenRTRBase);
   addpath(GenRTRBase);
   addpath([GenRTRBase '/drivers']);
   addpath([GenRTRBase '/solvers']);
   addpath([GenRTRBase '/tests']);
