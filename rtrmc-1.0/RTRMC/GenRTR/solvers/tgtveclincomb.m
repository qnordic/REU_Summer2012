function w = tgtveclincomb(a1, v1, a2, v2)

    if nargin == 2
    
        if ~iscell(v1)
            w = a1*v1;
        else
            w = cell(size(v1));
            for i = 1 : length(w)
                w{i} = a1*v1{i};
            end
        end
        
    elseif nargin == 4

        if ~iscell(v1)
            w = a1*v1 + a2*v2;
        else
            w = cell(size(v1));
            for i = 1 : length(w)
                w{i} = a1*v1{i} + a2*v2{i};
            end
        end
        
    else
        
        error('Incorrect use of arguments');
        
    end
    
end