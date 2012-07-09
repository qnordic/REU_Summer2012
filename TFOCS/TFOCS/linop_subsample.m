function op = linop_subsample( sz, omega, SYMMETRIC )
%OP = LINOP_SUBSAMPLE( SZ, OMEGA )
%   vector and matrix subsampling. Depending on SZ and OMEGA,
%   this can do row-sampling (e.g. a partial FFT)
%   or it can sample specific entries of a matrix (e.g. matrix completion)
%
%OP = LINOP_SUBSAMPLE( OMEGA )
%   works if OMEGA is a sparse matrix, whos nonzero entries
%   specify the entries to sample.
%
%OP = LINOP_SUBSAMPLE( ..., SYMMETRIC )
%   If SYMMETRIC is true, then
%   forces the domain to be the space of symmetric matrices


% Designed to be used with a fft or dct
% Do we need a separate oversampling version?
% Also, make linop_DCT and linop_FFT ?
%   The reason we might want this is that for linop_FFT,
%   people will probably use idct as the transpose -- this is not 
%   correct, due to scaling

% Help documentation: TBD
%    Constructs a TFOCS-compatible linear operator from separate function
%    handles that compute the forward and adjoint operations. The first
%    argument, SZ, gives the size of the linear operator; and the forward
%    and adjoint handles are AF and AT, respectively.
% 
%    If the inputs and outputs are simple vectors, then SZ can take the
%    standard Matlab form [N,M], where N is the input length and M is the
%    output length. If the input or output is a matrix or array, then SZ
%    should take the form { S_in, S_out }, where S_in is the size of the
%    input and S_out is the size of the output.
%

error(nargchk(1,3,nargin));
if nargin < 3
    SYMMETRIC = false;
end
if nargin == 1
    omega = sz;
    sz = { size(omega), [nnz(omega),1] };
elseif nargin ==2 && issparse(sz)
    SYMMETRIC = omega;
    omega = sz;
    sz = { size(omega), [nnz(omega),1] };
end

if numel( sz ) ~= 2,
    error( 'Size must have two elements.' );
elseif isnumeric( sz ),
    sz = { [sz(2),1], [sz(1),1] };
elseif ~isa( sz, 'cell' ),
    error( 'Invalid operator size specification.' );
end

dom = sz{1};
n1 = dom(1); n2 = dom(2);
ran = sz{2};

% There are several possibilities.  Let x be a point in the domain
%{
    x is a vector.  Then omega should be a vector.
                    We return x(omega), resp.

    x is a matrix.
        omega is a vector
            This is ambiguous: does the user want x(omega,:) or x(omega)?
        omega is a matrix with 2 columns, then assume it is [I,J]
            We convert it to linear indices.
        omega is a general matrix, more than 2 columns
            Not sure what the user means; report an error.
        omega is a sparse matrix
            We find it's entries, and use those

%}

if n2 == 1
    % x is a vector, not a matrix.  Simple.
    op = @(x,mode) linop_subsample_vector( sz, omega, x, mode );
else
    % trickier case.
    if issparse(omega)
        ind = find(omega);
        [I,J] = ind2sub( sz{1}, ind );
        op = @(x,mode) linop_subsample_matrix( sz, ind, I, J, SYMMETRIC, x, mode );
    elseif isvector(omega)
        op = @(x,mode) linop_subsample_matrix( sz, omega, SYMMETRIC,x, mode );
    elseif size(omega,2) == 2
        indx = sub2ind( sz{1}, omega(:,1), omega(:,2) );
        op = @(x,mode) linop_subsample_matrix( sz, indx, SYMMETRIC,x, mode );
    else
        error('omega is not an acceptable size; perhaps you meant it to be sparse?');
    end
        
    
end


function y = linop_subsample_vector(sz, omega, x, mode )
switch mode,
    case 0, y = sz;
    case 1, y = x(omega);
    case 2, 
        y = zeros( sz{1} );
        y(omega) = x;
end
function y = linop_subsample_matrix(sz, omega, indI, indJ,SYMMETRIC, x, mode )
switch mode,
    case 0, y = sz;
    case 1
        S = [];
        S.type = '()';
        S.subs = {omega};
        y = subsref(x,S);
    case 2, 
        dom = sz{1}; n1 = dom(1); n2 = dom(2);
        y = sparse( indI, indJ, x, n1, n2 );
        if SYMMETRIC
            % in future, might update this
            % e.g. force omega to only refer to lower part of matrx,
            % and then do the update y = y + tril(y,-1)'
            y = (y+y')/2;
        end
end



% TFOCS v1.1 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
