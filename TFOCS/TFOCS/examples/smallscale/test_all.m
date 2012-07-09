%{
Runs all tests in this directory.

Recommended if you want to be sure that TFOCS is working 100% perfectly
on your system.

The scripts will produce an error message if anything goes wrong.
So if you don't seen any error message in red text, then everything
is working.

%}

w   = what;
fileList    = w.m;
% and exclude the current file, otherwise we go into an infinite loop:
fileList    = fileList(  ~strcmpi( fileList, 'test_all.m' ) );
tm_all = tic;
for k = 1:length(fileList)
    
    file    = fileList{k}; % this has the .m extension, which we need to remove
    eval( file(1:end-2) );
    
end
fprintf('\n\n\n-- Success!! --\n');

toc( tm_all )

% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
