function flat = flattenDynamicMat(M,n,symb)
    assert(all(size(M)==[6,6,4,4]))
    if nargin==2 || symb
        flat = sym(zeros(7*n));
    else 
        flat = zeros(7*n);
    end
    for i = 1:n
        for j = 1:n
            flat(7*i-6,7*j-6) = M(1,1,2*i-1,2*j-1);
            flat(7*i-5:7*i,7*j-5:7*j) = M(:,:,2*i,2*j);
            flat(7*i-6,7*j-5:7*j) = M(1,:,2*i-1,2*j);
            flat(7*i-5:7*i,7*j-6) = M(:,1,2*i,2*j-1);
        end
    end
end