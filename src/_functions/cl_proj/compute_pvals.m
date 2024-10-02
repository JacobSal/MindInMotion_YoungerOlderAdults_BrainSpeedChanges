function pvals = compute_pvals(oridat, surrog, tail)
    
    if nargin < 3
        tail = 'both';
    end;
    
    if myndims(oridat) > 1        
        if size(oridat,2) ~= size(surrog, 2) | myndims(surrog) == 2
            if size(oridat,1) == size(surrog, 1)
                surrog = repmat( reshape(surrog, [size(surrog,1) 1 size(surrog,2)]), [1 size(oridat,2) 1]);
            elseif size(oridat,2) == size(surrog, 1)
                surrog = repmat( reshape(surrog, [1 size(surrog,1) size(surrog,2)]), [size(oridat,1) 1 1]);
            else
                error('Permutation statistics array size error');
            end;
        end;
    end;

    surrog = sort(surrog, myndims(surrog)); % sort last dimension
    
    if myndims(surrog) == 1    
        surrog(end+1) = oridat;        
    elseif myndims(surrog) == 2
        surrog(:,end+1) = oridat;        
    elseif myndims(surrog) == 3
        surrog(:,:,end+1) = oridat;
    else
        surrog(:,:,:,end+1) = oridat;
    end;

    [tmp idx] = sort( surrog, myndims(surrog) );
    [tmp mx]  = max( idx,[], myndims(surrog));        
                
    len = size(surrog,  myndims(surrog) );
    pvals = 1-(mx-0.5)/len;
    if strcmpi(tail, 'both')
        pvals = min(pvals, 1-pvals);
        pvals = 2*pvals;
    end;    
    
    function val = myndims(a)
    if ndims(a) > 2
        val = ndims(a);
    else
        if size(a,1) == 1,
            val = 2;
        elseif size(a,2) == 1,
            val = 1;
        else
            val = 2;
        end;
    end; 
