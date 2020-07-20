function flat = flattenStateMat(q,n,symb)
    if nargin==2 || symb
        flat = sym(zeros(7*n,1));
    else 
        flat = zeros(7*n,1);
    end
    for i = 1:n
        %the 1x1s
        flat(7*i-6) = q(2*i-1,1);
        %the 6x6s
        flat(7*i-5:7*i) = q(2*i,:);
    end
end