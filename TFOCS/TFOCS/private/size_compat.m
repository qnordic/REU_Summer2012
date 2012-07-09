function [ a, sX ] = size_compat( sX, sY )
a = true;
switch class( sX ),
    case 'double',
        if isempty( sX ) || all( sX == 0 ),
            sX = sY;
        elseif isempty( sY ) || all( sY == 0 ),
        elseif ~isequal( sX, sY ),
            
            % Feb 29, 2012. Special case:
            %   One reprsents the size a x b x c, where c = 1
            %   The other is a x b (since Matlab often automatically squeezes
            %   3D arrays to 2D if the 3rd dimension is a singletone)
            if (length(sX) >= 3 && length(sX) == length(sY)+1 && sX(end)==1) || ...
                    (length(sY) >= 3 && length(sY) == length(sX)+1 && sY(end)==1)
                % do nothing
            else
                a = false;
            end
        end
    case 'cell',
        if ~isa( sY, 'cell' ) || numel( sX ) ~= numel( sY ) || isa( sX{1}, 'function_handle' ) && ~isequal( sX, sY ),
            a = false;
        elseif isa( sX{1}, 'function_handle' ),
            a = isequal( sX, sY );
        else
            for k = 1 : numel( sX ),
                [ ta, sX{k} ] = size_compat( sX{k}, sY{k} );
                a = a && ta;
            end
        end
    otherwise,
        a = isequal( sX, sY );
end
if ~a,
    sX = [];
end

% TFOCS v1.1a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

