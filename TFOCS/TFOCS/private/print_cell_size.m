function print_cell_size( c , fid, offsetPossible )
% c should be a cell array

if nargin < 2 || isempty(fid), fid = 1; end
if nargin < 3 || isempty(offsetPossible), offsetPossible = false; end

if ~iscell(c)
    fprintf(fid,'Input to printCellSize should be a cell array');
    return;
else
    for j = 1:length(c)
        if j == length(c) && offsetPossible && all(size( c{j} ) == [1,2] ) ...
                && all( c{j} == [1,1] )
            fprintf(fid,'\tcomponent %2d is fixed (i.e. an offset)\n', j );
        else
            fprintf(fid,'\tcomponent %2d: ', j );
            if isempty( c{j} )
                fprintf('size not yet determined\n');
            else
                d = c{j};
                if length(d) < 2, d = [d,1]; end % this case shouldn't arise...
                for k = 1:(length(d)-1)
                    fprintf('%4d x ',d(k) );
                end
                fprintf('%4d\n', d(k+1) ); % bug, Feb 29 2012: change d(k) to d(k+1)
            end
        end
    end
end

% TFOCS v1.1a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.