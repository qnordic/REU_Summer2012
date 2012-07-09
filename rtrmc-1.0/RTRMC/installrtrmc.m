% RTRMC installation script.
% 
% BEFORE your run this script, try to simply execute TestRTRMC and see if
% it works. If it doesn't, then it is probably because the C-Mex codes need
% to be compiled. Launching this script should do the trick. If compilation
% fails, then you probably need to set up your C-compiler correctly for
% Matlab. Here are the official instructions for that:
% http://www.mathworks.nl/support/compilers/R2011b/win32.html
% From there, instructions for different Matlab versions and operating
% systems are easily reachable.
%
% if you have trouble installing/using this code, feel free to contact the
% author at: nicolas.boumal@uclouvain.be
%
% Nicolas Boumal, UCLouvain, May 19, 2011.

% Put this flag to 'true' if TestRTRMC failed.
I_launched_TestRTRMC_and_it_failed = true;

if ~I_launched_TestRTRMC_and_it_failed
    error(['Please first try to launch the script TestRTRMC. ' ...
       'If TestRTRMC executes without errors, then '...
       'there is no need to launch installrtrmc. '...
       'If it fails, then edit the flag on line 15 and launch the '...
       'present script again. ']);
end
   
mex -largeArrayDims spmaskmult.c
mex -largeArrayDims setsparseentries.c

if ~isempty(strfind(computer, 'WIN'))
    mex -largeArrayDims -lmwlapack cholsolvecell.c
    mex -largeArrayDims -lmwlapack -lmwblas spbuildmatrix.c
    mex -largeArrayDims -lmwlapack -lmwblas buildAchol.c
else
    mex -largeArrayDims -llapack cholsolvecell.c
    mex -largeArrayDims -llapack -lblas spbuildmatrix.c
    mex -largeArrayDims -llapack -lblas buildAchol.c
end
