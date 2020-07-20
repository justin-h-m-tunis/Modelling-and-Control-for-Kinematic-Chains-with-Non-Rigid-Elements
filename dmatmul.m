function M = dmatmul(A,B, sym, dimA, dimB)
    na = ndims(A);
    nb = ndims(B);
    assert(nb >= na, 'B must have the same or greater number of dims than A')
    if nargin==5
        perma = [dimA(1),dimA(2), find((1:na ~= dimA(1)) & (1:na ~= dimA(2)))];
        permb = [dimB(1),dimB(2), find((1:nb ~= dimB(1)) & (1:nb ~= dimB(2)))];
        A = permute(A,perma);
        B = permute(B,permb);
        na = ndims(A);
        nb = ndims(B);
    end
    outshape = size(B);
    outshape(1) = size(A,1);
    M = zeros(outshape);
    if nargin==2
        sym = 1;
    end
    if sym == 1
        if max(na, nb) == 2
            M = A * B;
            return
        else
            subsA.type='()';
            for i = 1:na-1
                subsA.subs{i} = ':';
            end
            subsB.type='()';
            for i = 1:nb-1
                subsB.subs{i} = ':';
            end
            last_dim = size(B,nb);
            for i = 1:last_dim
                subsB.subs{nb} = i;
                if na==nb && size(A,nb) ~= 1
                   subsA.subs{na} = i ;
                else
                   subsA.subs{na} = ':';
                end
                M = subsasgn(M,subsB,dmatmul(subsref(A,subsA),subsref(B,subsB)));
            end
        end
    else
        pagefun(@mtimes,A,B);
    end
end
